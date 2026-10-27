%% metabolomicsOAT.m
% Thesis Ch5: which personalisation inputs (blood metabolites, physiological
% parameters) change which model predictions, on one common scale?
%
% Each input is changed ON ITS OWN between the low and high end of its healthy range,
% everything else at the control (the fed Harvey/Harvetta of fastingWBM.m, also the
% control of metabolomicsPersonalisation.m):
%   metabolites  2.5th / 97.5th percentile of the healthy distribution per sex
%                (healthyReferencePanel.m), personalised with the persWBMmetabolomics
%                window (+/-20%) through metabBoundMap / applyMetabBounds (verified
%                identical to persWBMmetabolomics in metabolomicsPersonalisation.m Part 4)
%   physiology   healthy range per sex (healthyPhysiologyRanges.m), rebuilt with the
%                fastingWBM.m pipeline; body weight also recomputes organ fractions
% Predictions = FVA min/max of 11 organ exchange read-outs under BMR conditions
% (fastingTargetFVA). LP optimal values are unique, so the metric is free of the
% alternative-optima noise that affects distances between flux vectors.
%   effect(input, read-out) = (|max_hi - max_lo| + |min_hi - min_lo|) / control range
%
% Parallel sessions (start them a few minutes apart: initCobraToolbox updates git
% submodules and two sessions starting together collide):
%   sessionId = 1; nSessions = 3; metabolomicsOAT
% Each job saves its own file in results/metabolomics/oatFVA; rerun to resume. Session 1
% also runs two checks: the control rebuild, and whether the plain min-norm QP (no lifting)
% solves on the control. Run with nSessions = 0 to only combine finished jobs into the
% ranking, heatmap and tables.
%
% Anna Sheehy 2026.

% oatMode = 'l2' instead solves each case once for its basal flux distribution (min-norm QP,
% calculateModelBMR: lifted coupling + barrier) and saves BMR and the whole flux vector,
% in results/metabolomics/oatL2. Height changes no bounds, so its low/high pair gives the
% numerical noise floor used to decide which reactions changed.

%% Setup ------------------------------------------------------------------------------
if ~exist('sessionId', 'var'), sessionId = 1; end
if ~exist('nSessions', 'var'), nSessions = 1; end
if ~exist('oatMode', 'var'), oatMode = 'fva'; end
win    = 0.2;
% Physiological inputs are set to the low and high ends of their healthy range
% (references/healthyPhysiologyRanges.csv, healthyPhysiologyRanges.m), as the metabolites
% are set to the 2.5th and 97.5th percentiles of theirs.
% HeightBSA: height acting through cardiac output = cardiac index x body surface area
% (Du Bois & Du Bois 1916), with the cardiac index fixed so the reference individual
% keeps its cardiac output (by default height changes only blood volume, which no
% constraint uses).
physInputs = {'bodyWeight', 'Body weight'; 'Height', 'Height'; 'HeartRate', 'Heart rate'; ...
              'Hematocrit', 'Haematocrit'; 'GlomerularFiltrationRate', 'GFR'; ...
              'HeightBSA', 'Height (surface-area cardiac output)'; 'UrineCreatinine', 'Urine creatinine'};
readout = oatReadouts();                                          % 11 organ exchanges {rxn, label}

P = persWBMpaths();
fastPath  = fullfile(P.results, 'fasting');
resPath   = fullfile(P.results, 'metabolomics');
outDir    = fullfile(resPath, ['oat' upper(oatMode)]);
figPath   = fullfile(resPath, 'figures');
thesisFig = P.figures;
panelFile = fullfile(P.references, 'healthyReferenceValues.csv');
physFile  = fullfile(P.references, 'healthyPhysiologyRanges.csv');
if ~exist(outDir, 'dir'), mkdir(outDir); end
if ~exist(figPath, 'dir'), mkdir(figPath); end

if exist('startDelay', 'var'), pause(startDelay); end   % stagger initCobraToolbox across sessions
global CBT_LP_SOLVER
if isempty(CBT_LP_SOLVER), initCobraToolbox(false); end
changeCobraSolver('gurobi', 'LP', 0);
changeCobraSolver('gurobi', 'QP', 0);
if ~usejava('desktop'), set(0, 'DefaultFigureVisible', 'off'); end

sexes = {'male', 'female'};
S  = load(fullfile(fastPath, 'ckpt_feasibleModels.mat'), 'maleFed', 'femaleFed');
IPd = load(fullfile(fastPath, 'ckpt_HMDBconstrained.mat'), 'ipMale', 'ipFemale');
Rw = load(fullfile(fastPath, 'ckpt_rawModels.mat'), 'harvey', 'harvetta');
Mp = load(fullfile(resPath, 'ckpt_map.mat'), 'map');
ctrl = struct('male', S.maleFed, 'female', S.femaleFed);
ipar = struct('male', IPd.ipMale, 'female', IPd.ipFemale);
raw  = struct('male', Rw.harvey, 'female', Rw.harvetta);
map  = Mp.map;
clear S IPd Rw Mp
EUAverageDietNew;                                                 % defines Diet

panel  = healthyReferencePanel(panelFile);
panel  = panel(ismember(strcat(panel.vmh, '[bc]'), ctrl.male.mets) & ...
               ismember(strcat(panel.vmh, '[bc]'), ctrl.female.mets), :);
metIDs = panel.vmh(:)';  nMet = numel(metIDs);
% low/high levels per sex: 2.5th and 97.5th percentile of the healthy distribution (2 x nMet, uM)
levels = struct('male',   [panel.p025Male'; panel.p975Male'], ...
                'female', [panel.p025Female'; panel.p975Female']);
[physLim, physTab] = healthyPhysiologyRanges(physFile, ipar);   % physLim.(sex).(field) = [low high]
writetable(physTab, fullfile(resPath, 'physiologyRanges.csv'));
physiologyRangesTable(physFile, ipar, fullfile(P.tables, 'physiologyRanges.tex'));

% job list: {sex, kind, id, level}
jobs = cell(0, 4);
for k = 1:2
    jobs(end+1, :) = {sexes{k}, 'control', 0, 0}; %#ok<SAGROW>
    for q = 1:nMet, for L = 1:2, jobs(end+1, :) = {sexes{k}, 'metabolite', q, L}; end, end %#ok<SAGROW>
    for p = 1:size(physInputs, 1), for L = 1:2, jobs(end+1, :) = {sexes{k}, 'physiology', p, L}; end, end %#ok<SAGROW>
end
jobFile = @(j) fullfile(outDir, sprintf('%s_%s_%02d_%d.mat', jobs{j, :}));

%% Checks (session 1) -------------------------------------------------------------------
chkFile = fullfile(outDir, 'checks.mat');
if strcmp(oatMode, 'fva') && nSessions > 0 && sessionId == 1 && ~exist(chkFile, 'file')
    chk = struct();
    for k = 1:2
        s = sexes{k};
        m = buildFedWBM(raw.(s), ipar.(s), Diet);                    % reference rebuild
        chk.(s).rebuildMaxDiff = max([abs(m.lb - ctrl.(s).lb); abs(m.ub - ctrl.(s).ub)]);
        % plain min-norm QP (no lifting) on the control
        A = ctrl.(s);  ob = strcmp(A.rxns, 'Whole_body_objective_rxn');
        A.c(:) = 0;  A.c(ob) = 1;  A.osenseStr = 'max';
        t0 = tic;  solLP = optimizeWBModel(A);
        m2 = A;  m2.lb(ob) = solLP.f * (1 - 1e-4);  m2.c(:) = 0;  m2.osenseStr = 'min';
        solQP = optimizeCbModel(m2, 'min', 1e-6, 1, struct('ScaleFlag', 0, 'feasTol', 1e-4));
        chk.(s).lpStat = solLP.stat;  chk.(s).qpStat = solQP.stat;
        chk.(s).qpOrigStat = solQP.origStat;  chk.(s).seconds = toc(t0);
        fprintf('Check %s: rebuild max |diff| %.3g; LP stat %d, plain min-norm QP stat %d (orig %s), %.0f s\n', ...
            s, chk.(s).rebuildMaxDiff, solLP.stat, solQP.stat, num2str(solQP.origStat), chk.(s).seconds);
    end
    save(chkFile, 'chk');
end

%% Jobs -------------------------------------------------------------------------------
if nSessions > 0
    myJobs = find(mod((1:size(jobs, 1)) - 1, nSessions) == sessionId - 1);
    if exist('reverseOrder', 'var') && reverseOrder, myJobs = fliplr(myJobs); end   % extra helper sessions
    if exist('onlyKinds', 'var'), myJobs = myJobs(ismember(jobs(myJobs, 2), onlyKinds)); end   % e.g. {'control','physiology'}
    for j = myJobs
        if exist(jobFile(j), 'file'), continue, end
        [s, kind, id, L] = jobs{j, :};
        t0 = tic;
        switch kind
            case 'control'
                m = ctrl.(s);  val = NaN;
            case 'metabolite'
                val = levels.(s)(L, id);
                m = applyMetabBounds(ctrl.(s), map.(s), metIDs(id), val, win);
            case 'physiology'
                fld = physInputs{id, 1};
                val = physLim.(s).(fld)(L);
                if strcmp(fld, 'HeightBSA')
                    ip0 = ipar.(s);
                    ip1 = recalcDependentParams(ip0, 'Height', val);
                    cardiacIndex = ip0.CardiacOutput / bsaDuBois(ip0.bodyWeight, ip0.Height);
                    ip1.CardiacOutput = cardiacIndex * bsaDuBois(ip1.bodyWeight, val);
                else
                    ip1 = recalcDependentParams(ipar.(s), fld, val);
                end
                orgW = 'none';  if strcmp(fld, 'bodyWeight'), orgW = 'std'; end   % body weight also sets organ weights
                m = onLiftedModel(ctrl.(s), buildFedWBM(raw.(s), ip1, Diet, orgW));
        end
        conflict = any(m.lb > m.ub + 1e-9);
        if strcmp(oatMode, 'fva')
            if conflict
                F = table(string(readout(:, 1)), nan(size(readout, 1), 1), nan(size(readout, 1), 1), ...
                    'VariableNames', {'rxn', 'fvaMin', 'fvaMax'});
            else
                F = fastingTargetFVA(m, readout(:, 1));
            end
            secs = toc(t0);
            save(jobFile(j), 'F', 'val', 'conflict', 's', 'kind', 'id', 'L', 'secs');
        else
            kcal = NaN;  stat = NaN;  v = [];
            if ~conflict
                b = calculateModelBMR(m);
                kcal = b.kcal;  stat = b.stat;  v = single(b.v);
            end
            secs = toc(t0);
            save(jobFile(j), 'kcal', 'stat', 'v', 'val', 'conflict', 's', 'kind', 'id', 'L', 'secs');
        end
        fprintf('Job %3d/%d %s %s %d L%d done (%.0f s)\n', j, size(jobs, 1), s, kind, id, L, secs);
    end
end

%% Combine ----------------------------------------------------------------------------
done = arrayfun(@(j) exist(jobFile(j), 'file') == 2, 1:size(jobs, 1));
if ~all(done)
    fprintf('Combine: %d/%d jobs finished; rerun with nSessions = 0 when all are done.\n', sum(done), numel(done));
    return
end
if strcmp(oatMode, 'l2'), combineL2(outDir, resPath, figPath, thesisFig, panel, physInputs, ctrl, sexes); return, end
nR = size(readout, 1);
inputs = [panel.name(:); physInputs(:, 2)];
kinds  = [repmat({'metabolite'}, nMet, 1); repmat({'physiology'}, size(physInputs, 1), 1)];
nIn = numel(inputs);
E = nan(nIn, nR, 2);  Dmax = nan(nIn, nR, 2);  Dmin = nan(nIn, nR, 2);  infeas = false(nIn, 2);
ctrlRange = nan(2, nR);
for k = 1:2
    s = sexes{k};
    C = load(fullfile(outDir, sprintf('%s_control_00_0.mat', s)));
    ctrlRange(k, :) = (C.F.fvaMax - C.F.fvaMin)';
    ctrlRange(k, ctrlRange(k, :) < 1e-9) = NaN;          % fixed read-out: effect undefined
    for i = 1:nIn
        if i <= nMet, kind = 'metabolite'; id = i; else, kind = 'physiology'; id = i - nMet; end
        lo = load(fullfile(outDir, sprintf('%s_%s_%02d_1.mat', s, kind, id)));
        hi = load(fullfile(outDir, sprintf('%s_%s_%02d_2.mat', s, kind, id)));
        Dmax(i, :, k) = (hi.F.fvaMax - lo.F.fvaMax)';
        Dmin(i, :, k) = (hi.F.fvaMin - lo.F.fvaMin)';
        E(i, :, k) = (abs(Dmax(i, :, k)) + abs(Dmin(i, :, k))) ./ ctrlRange(k, :);
        infeas(i, k) = any(isnan([lo.F.fvaMax; hi.F.fvaMax]));
    end
end
Emed = squeeze(median(E, 2, 'omitnan'));  Emax = squeeze(max(E, [], 2));
T = table(inputs, kinds, Emed(:, 1), Emed(:, 2), Emax(:, 1), Emax(:, 2), infeas(:, 1), infeas(:, 2), ...
    'VariableNames', {'input','kind','medianEffectMale','medianEffectFemale','maxEffectMale', ...
    'maxEffectFemale','infeasibleLevelMale','infeasibleLevelFemale'});
T.score = mean([T.maxEffectMale, T.maxEffectFemale], 2, 'omitnan');
[T, ord] = sortrows(T, 'score', 'descend');
writetable(T, fullfile(resPath, 'oatRanking.csv'));
long = [];
for k = 1:2
    [ii, rr] = ndgrid(1:nIn, 1:nR);
    Tk = table(repmat(string(sexes{k}), numel(ii), 1), string(inputs(ii(:))), string(kinds(ii(:))), ...
        string(readout(rr(:), 2)), reshape(E(:, :, k), [], 1), reshape(Dmin(:, :, k), [], 1), ...
        reshape(Dmax(:, :, k), [], 1), 'VariableNames', {'sex','input','kind','readout','effect', ...
        'deltaMin','deltaMax'});
    if isempty(long), long = Tk; else, long = [long; Tk]; end %#ok<AGROW>
end
writetable(long, fullfile(resPath, 'oatEffects_long.csv'));

fid = fopen(fullfile(resPath, 'summary_oat.txt'), 'w');
say = @(varargin) cellfun(@(h) fprintf(h, varargin{:}), {1, fid}, 'UniformOutput', false);
if exist(chkFile, 'file')
    load(chkFile, 'chk');
    for k = 1:2
        say('Check %s: control rebuild max |diff| %.3g; plain min-norm QP status %d\n', sexes{k}, ...
            chk.(sexes{k}).rebuildMaxDiff, chk.(sexes{k}).qpStat);
    end
end
say('Control FVA range per read-out (male): %s\n', strjoin(compose('%s %.1f', string(readout(:, 2)), ctrlRange(1, :)'), '; '));
say('Ranking (max effect over read-outs, mean of sexes; effect = shift in feasible range / control range):\n');
for r = 1:height(T)
    [~, top] = max(mean(E(ord(r), :, :), 3, 'omitnan'));
    say('  %2d %-24s %-10s  max M %.3f F %.3f  median M %.3f F %.3f  largest on %s%s\n', r, T.input{r}, ...
        T.kind{r}, T.maxEffectMale(r), T.maxEffectFemale(r), T.medianEffectMale(r), T.medianEffectFemale(r), ...
        readout{top, 2}, repmat(' [an infeasible level]', 1, any([T.infeasibleLevelMale(r), T.infeasibleLevelFemale(r)])));
end
fclose(fid);

% Figure: inputs x read-outs, Harvey and Harvetta (small metabolite effects in one row)
oatEffectsFigure(long, T, figPath, thesisFig, {'Urine creatinine'});   % reported in the text: its urea effect needs its own scale
fprintf('OAT: ranking in %s\n', fullfile(resPath, 'oatRanking.csv'));


%% Local function ----------------------------------------------------------------------
function combineL2(outDir, resPath, figPath, thesisFig, panel, physInputs, ctrl, sexes)
% Rank inputs by how much their low-to-high change alters the basal flux distribution
% (number of reactions changed, size of the change, organs involved) and BMR.
% A reaction counts as changed if its flux moves by more than tau = max(1e-3, twice the
% largest change between the two height cases), height changing no bounds.
nMet = height(panel);  nP = size(physInputs, 1);  nIn = nMet + nP;
inputs = [panel.name(:); physInputs(:, 2)];
kinds  = [repmat({'metabolite'}, nMet, 1); repmat({'physiology'}, nP, 1)];
nChg = nan(nIn, 2);  dist = nan(nIn, 2);  dBMR = nan(nIn, 2);  bad = false(nIn, 2);
topOrg = repmat({""}, nIn, 2);  tau = nan(1, 2);  ctrlBMR = nan(1, 2);
for k = 1:2
    s = sexes{k};
    org = regexp(ctrl.(s).rxns, '^([A-Za-z0-9]+)_', 'tokens', 'once');
    org(cellfun(@isempty, org)) = {{'exchange'}};
    org = cellfun(@(c) c{1}, org, 'UniformOutput', false);
    C = load(fullfile(outDir, sprintf('%s_control_00_0.mat', s)));
    ctrlBMR(k) = C.kcal;
    V = cell(nIn, 2);  K = nan(nIn, 2);
    for i = 1:nIn
        if i <= nMet, kind = 'metabolite'; id = i; else, kind = 'physiology'; id = i - nMet; end
        for L = 1:2
            R = load(fullfile(outDir, sprintf('%s_%s_%02d_%d.mat', s, kind, id, L)));
            if ~isempty(R.v) && ismember(R.stat, [1 3]), V{i, L} = double(R.v); K(i, L) = R.kcal; end
        end
    end
    h = nMet + find(strcmp(physInputs(:, 1), 'Height'));
    tau(k) = max(1e-3, 2 * max(abs(V{h, 2} - V{h, 1})));
    for i = 1:nIn
        if isempty(V{i, 1}) || isempty(V{i, 2}), bad(i, k) = true; continue, end
        dv = V{i, 2} - V{i, 1};
        chg = abs(dv) > tau(k);
        nChg(i, k) = sum(chg);  dist(i, k) = norm(dv);  dBMR(i, k) = K(i, 2) - K(i, 1);
        if any(chg)
            [u, ~, j] = unique(org(chg));  cnt = accumarray(j, 1);
            [cs, o] = sort(cnt, 'descend');  n3 = min(3, numel(cs));
            topOrg{i, k} = strjoin(compose('%s %d', string(u(o(1:n3))), cs(1:n3)), ', ');
        end
    end
end
T = table(inputs, kinds, nChg(:, 1), nChg(:, 2), dist(:, 1), dist(:, 2), dBMR(:, 1), dBMR(:, 2), ...
    string(topOrg(:, 1)), string(topOrg(:, 2)), bad(:, 1), bad(:, 2), 'VariableNames', ...
    {'input','kind','nChangedMale','nChangedFemale','distMale','distFemale','dBMRMale','dBMRFemale', ...
    'topOrgansMale','topOrgansFemale','noSolutionMale','noSolutionFemale'});
T.score = mean([T.nChangedMale, T.nChangedFemale], 2, 'omitnan');
T = sortrows(T, 'score', 'descend');
writetable(T, fullfile(resPath, 'oatL2Ranking.csv'));

fid = fopen(fullfile(resPath, 'summary_oatL2.txt'), 'w');
say = @(varargin) cellfun(@(f) fprintf(f, varargin{:}), {1, fid}, 'UniformOutput', false);
say('Control BMR (L2): male %.0f, female %.0f kcal/day. Change threshold tau: male %.3g, female %.3g mmol/day\n', ...
    ctrlBMR(1), ctrlBMR(2), tau(1), tau(2));
say('%-24s %-10s %8s %8s %10s %10s  %s\n', 'input', 'kind', 'nChg M', 'nChg F', 'dBMR M', 'dBMR F', 'organs (male)');
for r = 1:height(T)
    say('%-24s %-10s %8.0f %8.0f %10.1f %10.1f  %s%s\n', T.input{r}, T.kind{r}, T.nChangedMale(r), ...
        T.nChangedFemale(r), T.dBMRMale(r), T.dBMRFemale(r), T.topOrgansMale(r), ...
        repmat('  [a level had no solution]', 1, T.noSolutionMale(r) || T.noSolutionFemale(r)));
end
fclose(fid);

oatL2Figure(T, figPath, thesisFig);
end
