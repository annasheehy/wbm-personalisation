%% metabolomicsPersonalisation.m
% Thesis Ch5, Results 2.2: what does personalising a WBM with blood metabolomics do?
%
% Engine: physiologicalConstraintsHMDBbased(..., 'direct', ..., 'bc') (PSCM toolbox,
% Thiele), exactly as called by persWBMmetabolomics (Persephone, Nap et al.), unchanged.
% Inputs: SYNTHETIC blood metabolomes, log-normal draws per sex from healthy-adult
% reference studies (Trabado et al. 2017, n = 800; Nordic Reference Interval Project,
% Rustad et al. 2004; Hordaland Health Study, Konstantinova et al. 2008), transcribed in
% references/healthyReferenceValues.csv and fitted by healthyReferencePanel.m. Every
% result is a model output; only the input concentrations are synthetic.
%
% Each personalised bound is linear in the measured concentration (organ uptake
% lb = k*cmax; kidney filtration [k*cmax, k*cmin]; CSF [k*cmin, k*cmax]), so the
% coefficients k are read from two engine calls and individuals are built by scaling.
% Part 4 checks this against direct engine calls and against persWBMmetabolomics.
%
% PARTS (checkpointed in results/metabolomics; set redo.<name> = true to rerun one)
%   0 setup            1 models + panel        2 synthetic cohort
%   3 bound map        4 verification          5 feasibility x allowance, crossed bounds
%   6 figures + tables
%
% Anna Sheehy 2026.

%% ==========================================================================
%  PART 0. SETUP
%  ==========================================================================
nPerSex     = 25;              % synthetic individuals per sex
rngSeed     = 1;
windows     = [0.1 0.2 0.5];   % +/- measured concentration allowance (0.1 infant WBMs, 0.2 Persephone)
refWindow   = 0.2;             % Persephone default; used for the uptake-limit ratios (Fig B)
if ~exist('stopAfter', 'var'), stopAfter = 6; end   % last part to run (settable before calling)
redo = struct('cohort',false,'map',false,'verify',false,'feas',false);

P = persWBMpaths();
fastPath  = fullfile(P.results, 'fasting');
resPath   = fullfile(P.results, 'metabolomics');
figPath   = fullfile(resPath, 'figures');
thesisFig = P.figures;
panelFile = fullfile(P.references, 'healthyReferenceValues.csv');
if ~exist(figPath, 'dir'), mkdir(figPath); end
ck = @(name) fullfile(resPath, ['ckpt_' name '.mat']);

global CBT_LP_SOLVER
if isempty(CBT_LP_SOLVER), initCobraToolbox(false); end
changeCobraSolver('gurobi', 'LP', 0);
changeCobraSolver('gurobi', 'QP', 0);
if ~usejava('desktop'), set(0, 'DefaultFigureVisible', 'off'); end
fprintf('Part 0 - results: %s\n', resPath);

%% ==========================================================================
%  PART 1. MODELS + PANEL
%  ==========================================================================
% Controls = the fed, HMDB-constrained Harvey/Harvetta of fastingWBM.m Part 2. Persephone's
% control model is the same model passed through the engine with no data, which leaves
% every bound unchanged.
S = load(fullfile(fastPath, 'ckpt_feasibleModels.mat'), 'maleFed', 'femaleFed');
IPd = load(fullfile(fastPath, 'ckpt_HMDBconstrained.mat'), 'ipMale', 'ipFemale');
sexes = {'male', 'female'};
ctrl = struct('male', S.maleFed, 'female', S.femaleFed);
ipar = struct('male', IPd.ipMale, 'female', IPd.ipFemale);
clear S IPd
for s = sexes
    assert(isfield(ipar.(s{1}), 'sex') && strcmpi(ipar.(s{1}).sex, s{1}), ...
        'IndividualParameters.sex missing or wrong for %s', s{1});
end

panel = healthyReferencePanel(panelFile);
inModel = ismember(strcat(panel.vmh, '[bc]'), ctrl.male.mets) & ...
          ismember(strcat(panel.vmh, '[bc]'), ctrl.female.mets);
if any(~inModel)
    fprintf('Part 1 - dropped (no [bc] metabolite in both models): %s\n', strjoin(panel.vmh(~inModel), ', '));
end
panel = panel(inModel, :);
metIDs = panel.vmh(:)';
nMet = numel(metIDs);
[~, mi] = ismember(strcat(metIDs, '[bc]'), ctrl.male.mets);
fprintf('Part 1 - panel of %d metabolites (source name -> model metName):\n', nMet);
for q = 1:nMet
    fprintf('   %-8s %-24s -> %s\n', metIDs{q}, panel.name{q}, ctrl.male.metNames{mi(q)});
end
% log-normal parameters per sex (log uM), 1 x nMet
logMu  = struct('male', panel.muMale', 'female', panel.muFemale');
logSig = struct('male', panel.sigmaMale', 'female', panel.sigmaFemale');
if stopAfter < 2, return, end

%% ==========================================================================
%  PART 2. SYNTHETIC COHORT
%  ==========================================================================
% Log-normal per metabolite and sex (Part 1); independent draws (no source gives
% covariances). Where a source does not separate men and women, both sexes are drawn
% from the same distribution.
if ~exist(ck('cohort'), 'file') || redo.cohort
    rng(rngSeed);
    conc = struct();
    for s = sexes
        conc.(s{1}) = exp(logMu.(s{1}) + logSig.(s{1}) .* randn(nPerSex, nMet));   % uM, nPerSex x nMet
    end
    save(ck('cohort'), 'conc', 'metIDs', 'rngSeed');
    % the same cohort in Persephone metadata format (ID, Sex, metabolites, units row)
    ids = cellstr([compose('M%02d', 1:nPerSex), compose('F%02d', 1:nPerSex)]');
    sexCol = [repmat({'male'}, nPerSex, 1); repmat({'female'}, nPerSex, 1)];
    body = [ids, sexCol, num2cell([conc.male; conc.female])];
    meta = [[{'ID', 'Sex'}, metIDs]; [{'unit', ''}, repmat({'umol/L'}, 1, nMet)]; body];
    writecell(meta, fullfile(resPath, 'syntheticCohort_persephoneFormat.csv'));
else
    load(ck('cohort'), 'conc');
end
fprintf('Part 2 - cohort: %d per sex, %d metabolites\n', nPerSex, nMet);
if stopAfter < 3, return, end

%% ==========================================================================
%  PART 3. BOUND MAP: which bounds one measured concentration sets, and how
%  ==========================================================================
if ~exist(ck('map'), 'file') || redo.map
    map = struct();
    for s = sexes
        t0 = tic;
        map.(s{1}) = metabBoundMap(ctrl.(s{1}), ipar.(s{1}), metIDs);
        fprintf('Part 3 - %s: %d bounds on %d reactions (%.0f s)\n', s{1}, height(map.(s{1})), ...
            numel(unique(map.(s{1}).rxnIdx)), toc(t0));
    end
    save(ck('map'), 'map');
else
    load(ck('map'), 'map');
end
if stopAfter < 4, return, end

%% ==========================================================================
%  PART 4. VERIFICATION
%  ==========================================================================
% (a) linear reconstruction == engine, for 3 individuals per sex and every window
% (b) engine == persWBMmetabolomics itself, with mixed input units (mmol/L, mg/dL, umol/L)
if ~exist(ck('verify'), 'file') || redo.verify
    ver = table('Size', [0 4], 'VariableTypes', {'string','double','double','double'}, ...
        'VariableNames', {'sex','individual','window','maxAbsDiff'});
    for s = sexes
        for i = 1:3
            for w = windows
                c = conc.(s{1})(i, :);
                dataBC = [{'met','minCon','maxCon'}; [metIDs(:), num2cell(c(:)*(1-w)), num2cell(c(:)*(1+w))]];
                mEng = physiologicalConstraintsHMDBbased(ctrl.(s{1}), ipar.(s{1}), '', 'direct', dataBC, 'bc');
                mLin = applyMetabBounds(ctrl.(s{1}), map.(s{1}), metIDs, c, w);
                d = max([abs(mEng.lb - mLin.lb); abs(mEng.ub - mLin.ub)]);
                ver(end+1, :) = {string(s{1}), i, w, d}; %#ok<SAGROW>
            end
        end
    end
    % (b) persWBMmetabolomics with mixed units for glucose, creatinine, alanine
    persCheck = struct('ok', false, 'maxRelDiff', NaN, 'message', "");
    try
        demo = {'glc_D', 'crtn', 'ala_L'};
        c = conc.male(1, :);
        cv = c(cellfun(@(x) find(strcmp(metIDs, x)), demo));
        MWcrtn = getMolecularMass('C4H7N3O');               % creatinine
        vals = {cv(1) / 1000, cv(2) * MWcrtn / 1e4, cv(3)};  % mmol/L, mg/dL, umol/L
        T = cell2table([{'unit'; 'value'; 'compartment'}, ...
            [{'mmol/L', 'mg/dL', 'umol/L'}; vals; {'blood', 'blood', 'blood'}]], ...
            'VariableNames', [{'ID'}, demo]);
        mP = ctrl.male;  mP.ID = 'verify_M01';  mP.SetupInfo.IndividualParameters = ipar.male;
        mP = persWBMmetabolomics('male', T, 'iWBM', mP);
        dataBC = [{'met','minCon','maxCon'}; [demo(:), num2cell(cv(:)*0.8), num2cell(cv(:)*1.2)]];
        mEng = physiologicalConstraintsHMDBbased(ctrl.male, ipar.male, '', 'direct', dataBC, 'bc');
        chg = find(mEng.lb ~= ctrl.male.lb | mEng.ub ~= ctrl.male.ub);
        rel = max(abs([mP.lb(chg) - mEng.lb(chg); mP.ub(chg) - mEng.ub(chg)]) ./ ...
                  max(abs([mEng.lb(chg); mEng.ub(chg)]), 1e-12));
        persCheck = struct('ok', rel < 1e-3, 'maxRelDiff', rel, ...
            'message', sprintf('%d bounds compared', numel(chg)));
    catch err
        persCheck.message = string(err.message);
    end
    save(ck('verify'), 'ver', 'persCheck');
else
    load(ck('verify'), 'ver', 'persCheck');
end
fprintf('Part 4 - linear vs engine: max abs diff %.3g over %d checks\n', max(ver.maxAbsDiff), height(ver));
fprintf('Part 4 - persWBMmetabolomics (mixed units) vs engine: ok=%d, max rel diff %.3g (%s)\n', ...
    persCheck.ok, persCheck.maxRelDiff, persCheck.message);
if stopAfter < 5, return, end

%% ==========================================================================
%  PART 5. FEASIBILITY x WINDOW (+ culprit metabolites)
%  ==========================================================================
% Persephone criterion: unit whole-body maintenance flux achievable. Also tested with
% the resting brain/heart ATP floors used for BMR and the read-outs.
if ~exist(ck('feas'), 'file') || redo.feas
    feas = table('Size', [0 7], 'VariableTypes', {'string','double','double','double','logical','logical','double'}, ...
        'VariableNames', {'sex','individual','window','nConflict','feasible','feasibleFloors','seconds'});
    ctrlFeas = struct();
    for s = sexes
        ctrlFeas.(s{1}) = [isFeasibleWBM(ctrl.(s{1}), false), isFeasibleWBM(ctrl.(s{1}), true)];
        for w = windows
            for i = 1:nPerSex
                t0 = tic;
                m = applyMetabBounds(ctrl.(s{1}), map.(s{1}), metIDs, conc.(s{1})(i, :), w);
                nC = sum(m.lb > m.ub + 1e-9);
                f1 = nC == 0 && isFeasibleWBM(m, false);
                f2 = f1 && isFeasibleWBM(m, true);
                feas(end+1, :) = {string(s{1}), i, w, nC, f1, f2, toc(t0)}; %#ok<SAGROW>
                fprintf('Part 5 - %s %02d w=%.1f conflicts=%d feasible=%d floors=%d (%.0f s)\n', ...
                    s{1}, i, w, nC, f1, f2, toc(t0));
            end
        end
    end
    save(ck('feas'), 'feas', 'ctrlFeas');
else
    load(ck('feas'), 'feas', 'ctrlFeas');
end
fprintf('Part 5 - median LP check time %.1f s per model\n', median(feas.seconds));

%% --------------------------------------------------------------------------
%  PART 5b. WHY MODELS FAIL, AND A FIX
%  --------------------------------------------------------------------------
% For every infeasible model: (a) which reactions have lb > ub and which personalised
% metabolite set them; (b) without a conflict, which single metabolite's restoration to
% its default bounds restores feasibility (culprit); (c) fix test, "measurement takes
% precedence": where a personalised bound crosses a fixed bound, the fixed bound is
% relaxed only as far as needed (moved onto the personalised one), then feasibility is
% re-tested. No diet or other compensation is applied anywhere.
if ~exist(ck('conflicts'), 'file') || redo.feas
    conflicts = table('Size', [0 7], 'VariableTypes', {'string','double','double','string','string','double','double'}, ...
        'VariableNames', {'sex','individual','window','rxn','met','personalised','fixedBound'});
    culprit = table('Size', [0 4], 'VariableTypes', {'string','double','double','string'}, ...
        'VariableNames', {'sex','individual','window','met'});
    fixTest = table('Size', [0 5], 'VariableTypes', {'string','double','double','logical','logical'}, ...
        'VariableNames', {'sex','individual','window','feasible','feasibleFloors'});
    for s = sexes
        mp = map.(s{1});
        bad = feas(feas.sex == s{1} & ~feas.feasibleFloors, :);
        for r = 1:height(bad)
            i = bad.individual(r);  w = bad.window(r);
            m = applyMetabBounds(ctrl.(s{1}), mp, metIDs, conc.(s{1})(i, :), w);
            cf = find(m.lb > m.ub + 1e-9);
            if ~isempty(cf)
                for j = cf'
                    row = mp(mp.rxnIdx == j, :);
                    if row.bound(1) == "lb", pv = m.lb(j); fb = m.ub(j); else, pv = m.ub(j); fb = m.lb(j); end
                    conflicts(end+1, :) = {string(s{1}), i, w, string(m.rxns{j}), string(row.met{1}), pv, fb}; %#ok<SAGROW>
                end
                mFix = m;
                for j = cf'
                    row = mp(mp.rxnIdx == j, :);
                    if row.bound(1) == "lb", mFix.ub(j) = mFix.lb(j); else, mFix.lb(j) = mFix.ub(j); end
                end
                f1 = isFeasibleWBM(mFix, false);  f2 = f1 && isFeasibleWBM(mFix, true);
                fixTest(end+1, :) = {string(s{1}), i, w, f1, f2}; %#ok<SAGROW>
                fprintf('Part 5b - %s %02d w=%.1f: %d conflict(s); after fix feasible=%d floors=%d\n', ...
                    s{1}, i, w, numel(cf), f1, f2);
            else
                for q = 1:nMet
                    keep = ~strcmp(mp.met, metIDs{q});
                    m2 = applyMetabBounds(ctrl.(s{1}), mp(keep, :), metIDs, conc.(s{1})(i, :), w);
                    if isFeasibleWBM(m2, true)
                        culprit(end+1, :) = {string(s{1}), i, w, string(metIDs{q})}; %#ok<SAGROW>
                    end
                end
                fprintf('Part 5b - %s %02d w=%.1f: no conflict; culprits: %s\n', s{1}, i, w, ...
                    strjoin(culprit.met(culprit.sex == s{1} & culprit.individual == i & culprit.window == w), ', '));
            end
        end
    end
    save(ck('conflicts'), 'conflicts', 'culprit', 'fixTest');
else
    load(ck('conflicts'), 'conflicts', 'culprit', 'fixTest');
end
writetable(conflicts, fullfile(resPath, 'conflicts.csv'));
writetable(culprit,   fullfile(resPath, 'culprits.csv'));
writetable(fixTest,   fullfile(resPath, 'fixTest.csv'));

% Conflict thresholds, exact rather than sampled: each personalised bound that faces a
% fixed bound of the control crosses it at one concentration (bounds are linear in c);
% report that concentration and the share of the healthy distribution beyond it.
thresholds = table('Size', [0 7], 'VariableTypes', {'string','string','string','double','double','string','double'}, ...
    'VariableNames', {'sex','rxn','met','window','threshold_uM','conflictWhen','popFraction'});
Phi = @(x) 0.5 * erfc(-x / sqrt(2));
for s = sexes
    mp = map.(s{1});  c0 = ctrl.(s{1});
    bothPers = intersect(mp.rxnIdx(mp.bound == "lb"), mp.rxnIdx(mp.bound == "ub"));
    [~, qOf] = ismember(mp.met, metIDs);
    for r = find(~ismember(mp.rxnIdx, bothPers) & qOf > 0)'
        j = mp.rxnIdx(r);  isLb = mp.bound(r) == "lb";
        if isLb, fb = c0.ub(j); else, fb = c0.lb(j); end
        if abs(fb) >= 1e5, continue, end                     % effectively unbounded
        for w = windows
            a = mp.k(r) * ((1 + w) * strcmp(mp.dependsOn{r}, 'max') + (1 - w) * strcmp(mp.dependsOn{r}, 'min'));
            if a == 0 || fb / a <= 0, continue, end          % cannot cross for c > 0
            cStar = fb / a;
            below = xor(isLb, a > 0);                         % lb > ub (or ub < lb) for c < cStar?
            p = Phi((log(cStar) - logMu.(s{1})(qOf(r))) / logSig.(s{1})(qOf(r)));
            side = "below";
            if ~below, p = 1 - p;  side = "above"; end
            if p < 1e-4, continue, end
            thresholds(end+1, :) = {string(s{1}), string(mp.rxn{r}), string(mp.met{r}), w, cStar, side, p}; %#ok<SAGROW>
        end
    end
end
writetable(thresholds, fullfile(resPath, 'conflictThresholds.csv'));
for r = 1:height(thresholds)
    fprintf('Part 5b - threshold %s %s (%s) w=%.1f: conflict %s %.0f uM, %.1f%% of the healthy distribution\n', ...
        thresholds.sex(r), thresholds.rxn(r), thresholds.met(r), thresholds.window(r), ...
        thresholds.conflictWhen(r), thresholds.threshold_uM(r), 100 * thresholds.popFraction(r));
end
if stopAfter < 6, return, end

%% ==========================================================================
%  PART 6. FIGURES + TABLES
%  ==========================================================================
typeOrder = {'organ uptake', 'kidney filtration', 'blood-brain barrier'};
sumFile = fopen(fullfile(resPath, 'summary.txt'), 'w');
say = @(varargin) cellfun(@(f) fprintf(f, varargin{:}), {1, sumFile}, 'UniformOutput', false);

% --- Table: the healthy ranges behind the synthetic profiles ------------------------
% HMDB upper limit = concentration behind the control's default organ uptake bound
% (default bound = k * c, so c = default / k; median over organs).
cDef = nan(nMet, 1);
mp = map.male;
sel = strcmp(mp.type, 'organ uptake') & mp.bound == "lb" & mp.defaultValue < 0;
for q = 1:nMet
    r = sel & strcmp(mp.met, metIDs{q});
    cDef(q) = median(mp.defaultValue(r) ./ mp.k(r));
end
Tin = panel(:, {'name', 'vmh', 'source', 'sexSpecific', 'p025Male', 'p975Male', 'p025Female', 'p975Female'});
Tin.hmdbUpperLimit = cDef;
tablePath = P.tables;
if ~exist(tablePath, 'dir'), mkdir(tablePath); end
% flag: HMDB upper limit below the 97.5th percentile of the healthy distribution (either sex)
lowHMDB = cDef < max(panel.p975Male, panel.p975Female);
% metabolites without a usable HMDB value get the toolbox default blood concentration
isDefault = abs(cDef - ipar.male.MConDefaultBc) < 1e-9;
Tin.hmdbBelowHealthy = lowHMDB;  Tin.toolboxDefault = isDefault;
writetable(Tin, fullfile(resPath, 'inputPanel.csv'));
% two tables: sexes reported together / separately (Tables/metabolomicsPanel.tex)
metabolomicsPanelTable(fullfile(resPath, 'inputPanel.csv'), fullfile(tablePath, 'metabolomicsPanel.tex'), ...
    nPerSex, ipar.male.MConDefaultBc);
say('Table: HMDB upper limit below the healthy 97.5th percentile for %d metabolites: %s\n', sum(lowHMDB), ...
    strjoin(panel.name(lowHMDB), ', '));

% --- Fig A: bounds set per metabolite, by type (identical metabolites as one bar) ----
for k = 1:2
    s = sexes{k};  mp = map.(s);
    cnt = zeros(nMet, numel(typeOrder));
    for q = 1:nMet
        for t = 1:numel(typeOrder)
            cnt(q, t) = numel(unique(mp.rxnIdx(strcmp(mp.met, metIDs{q}) & strcmp(mp.type, typeOrder{t}))));
        end
    end
    say('A %s: %d reactions (%d bounds); per metabolite median %d, range %d-%d\n', s, ...
        numel(unique(mp.rxnIdx)), height(mp), median(sum(cnt, 2)), min(sum(cnt, 2)), max(sum(cnt, 2)));
end
grpA = metabBoundMapFigure(map, metIDs, panel.name, figPath, thesisFig);
for j = find(cellfun(@numel, grpA.members) > 1)'
    say('A bar "%s": %s\n', grpA.label{j}, strjoin(grpA.members{j}, ', '));
end

% --- Fig B: personalised bound / default bound (organ uptake, refWindow) -------------
ratioMed = nan(nMet, 2);  Rall = struct();
for k = 1:2
    s = sexes{k};  mp = map.(s);
    sel = strcmp(mp.type, 'organ uptake') & mp.bound == "lb" & mp.defaultValue < 0 & isfinite(mp.defaultValue);
    R = nan(nPerSex, nMet);
    for i = 1:nPerSex
        m = applyMetabBounds(ctrl.(s), mp, metIDs, conc.(s)(i, :), refWindow);
        r = m.lb(mp.rxnIdx(sel)) ./ mp.defaultValue(sel);
        for q = 1:nMet
            R(i, q) = median(r(strcmp(mp.met(sel), metIDs{q})));
        end
    end
    Rall.(s) = R;
    ratioMed(:, k) = median(R, 1, 'omitnan')';
    say('B %s: median ratio <1 (tightened) for %d/%d metabolites; >1 for %d\n', s, ...
        sum(ratioMed(:, k) < 1), nMet, sum(ratioMed(:, k) > 1));
end
metabBoundRatioFigure(Rall, panel.name, refWindow, figPath, thesisFig);
writetable(table(panel.name, metIDs', ratioMed(:,1), ratioMed(:,2), 'VariableNames', ...
    {'metabolite','vmh','medianRatioMale','medianRatioFemale'}), fullfile(resPath, 'boundRatio.csv'));

% --- Fig C: feasibility by window, as personalised and after conflict resolution ----
pctC = metabFeasibilityFigure(feas, fixTest, windows, figPath, thesisFig);
for k = 1:2
    for wI = 1:numel(windows)
        r = feas(feas.sex == sexes{k} & feas.window == windows(wI), :);
        say(['C %s w=%.1f: feasible %.0f%%, with ATP floors %.0f%%, any conflict %d models; ' ...
             'after conflict resolution %.0f%% (%.0f%%)\n'], sexes{k}, windows(wI), pctC.(sexes{k})(wI, 1:2), ...
            sum(r.nConflict > 0), pctC.(sexes{k})(wI, 3:4));
    end
end
writetable(feas, fullfile(resPath, 'feasibility.csv'));
if ~isempty(conflicts)
    [u, ~, j] = unique(conflicts.rxn + " (" + conflicts.met + ")");
    say('C conflicting reactions (count of models x window): %s\n', strjoin(compose('%s x%d', u, accumarray(j, 1)), ', '));
end
for r = 1:height(thresholds)
    say('C threshold %s %s (%s) w=%.1f: conflict %s %.0f uM = %.1f%% of the healthy distribution\n', ...
        thresholds.sex(r), thresholds.rxn(r), thresholds.met(r), thresholds.window(r), ...
        thresholds.conflictWhen(r), thresholds.threshold_uM(r), 100 * thresholds.popFraction(r));
end
if ~isempty(fixTest)
    say('C fix "measurement takes precedence": %d/%d conflicting models feasible after fix (%d with ATP floors)\n', ...
        sum(fixTest.feasible), height(fixTest), sum(fixTest.feasibleFloors));
end
if ~isempty(culprit)
    [u, ~, j] = unique(culprit.met);
    say('C culprits (non-conflict infeasible models rescued by restoring one metabolite): %s\n', ...
        strjoin(compose('%s x%d', u, accumarray(j, 1)), ', '));
end

fclose(sumFile);
fprintf('Part 6 - figures in %s and %s\n', figPath, thesisFig);


%% ==========================================================================
%  LOCAL FUNCTIONS
%  ==========================================================================
function ok = isFeasibleWBM(m, withFloors)
% Unit whole-body maintenance flux achievable (Persephone criterion), optionally with
% the resting brain/heart ATP floors of calculateModelBMR.
if withFloors, m = setBMRConditions(m); end
ob = strcmp(m.rxns, 'Whole_body_objective_rxn');
m.c(:) = 0;  m.c(ob) = 1;  m.osenseStr = 'max';
sol = optimizeWBModel(m);
ok = sol.stat == 1 && sol.f >= 1 - 1e-6;
end
