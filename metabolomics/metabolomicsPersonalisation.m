%% metabolomicsPersonalisation.m
% Thesis Ch5, Results 2.2: what does personalising a WBM with serum metabolomics do?
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
%   3 bound map        4 verification          5 feasibility x window
%   6 read-out FVA     7 BMR (L2)              8 figures + tables
%
% Anna Sheehy 2026.

%% ==========================================================================
%  PART 0. SETUP
%  ==========================================================================
nPerSex     = 25;              % synthetic individuals per sex
rngSeed     = 1;
windows     = [0.1 0.2 0.5];   % +/- measurement window (0.1 infant WBMs, 0.2 Persephone)
refWindow   = 0.2;             % Persephone default; used for read-outs and BMR
maxFVA      = nPerSex;         % individuals per sex given read-out FVA
maxBMR      = nPerSex;         % individuals per sex given L2 BMR
if ~exist('stopAfter', 'var'), stopAfter = 8; end   % last part to run (settable before calling)
if ~exist('skip', 'var'), skip = []; end            % parts to skip, e.g. skip = [6 7] to draw figures early
redo = struct('cohort',false,'map',false,'verify',false,'feas',false, ...
              'fva',false,'bmr',false);

% Read-outs: organ exchanges benchmarked in the fasting chapter (fastingLiteratureTargets.m)
readout = {'Liver_EX_glc_D(e)_[bc]',     'Liver glucose'
           'Brain_EX_glc_D(e)_[csf]',    'Brain glucose'
           'Muscle_EX_glc_D(e)_[bc]',    'Muscle glucose'
           'RBC_EX_glc_D(e)_[bc]',       'RBC glucose'
           'Kidney_EX_glc_D(e)_[bc]',    'Kidney glucose'
           'Adipocytes_EX_glyc(e)_[bc]', 'Adipose glycerol'
           'Liver_EX_bhb(e)_[bc]',       'Liver 3-hydroxybutyrate'
           'Brain_EX_bhb(e)_[csf]',      'Brain 3-hydroxybutyrate'
           'Muscle_EX_ala_L(e)_[bc]',    'Muscle alanine'
           'Liver_EX_ala_L(e)_[bc]',     'Liver alanine'
           'EX_urea[u]',                 'Urinary urea'};

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
%  PART 6. READ-OUT FVA (refWindow, models feasible with floors)
%  ==========================================================================
if ismember(6, skip)
    fprintf('Part 6 - skipped\n');
elseif ~exist(ck('fva'), 'file') || redo.fva
    nR = size(readout, 1);
    fva = struct();
    for s = sexes
        F = fastingTargetFVA(ctrl.(s{1}), readout(:, 1));
        fva.(s{1}).ctrlMin = F.fvaMin';  fva.(s{1}).ctrlMax = F.fvaMax';
        fva.(s{1}).min = nan(nPerSex, nR);  fva.(s{1}).max = nan(nPerSex, nR);
        ok = feas.individual(feas.sex == s{1} & feas.window == refWindow & feas.feasibleFloors);
        for i = ok(1:min(end, maxFVA))'
            m = applyMetabBounds(ctrl.(s{1}), map.(s{1}), metIDs, conc.(s{1})(i, :), refWindow);
            F = fastingTargetFVA(m, readout(:, 1));
            fva.(s{1}).min(i, :) = F.fvaMin';  fva.(s{1}).max(i, :) = F.fvaMax';
            save(ck('fva'), 'fva');                       % incremental
            fprintf('Part 6 - FVA %s %02d done\n', s{1}, i);
        end
    end
    save(ck('fva'), 'fva');
else
    load(ck('fva'), 'fva');
end
if stopAfter < 7, return, end

%% ==========================================================================
%  PART 7. BMR (L2) at refWindow
%  ==========================================================================
if ismember(7, skip)
    fprintf('Part 7 - skipped\n');
elseif ~exist(ck('bmr'), 'file') || redo.bmr
    bmr = struct();
    for s = sexes
        b = calculateModelBMR(ctrl.(s{1}), 'norms', {'L2'});
        bmr.(s{1}).ctrl = b.L2.kcal;
        bmr.(s{1}).kcal = nan(nPerSex, 1);
        ok = feas.individual(feas.sex == s{1} & feas.window == refWindow & feas.feasibleFloors);
        for i = ok(1:min(end, maxBMR))'
            m = applyMetabBounds(ctrl.(s{1}), map.(s{1}), metIDs, conc.(s{1})(i, :), refWindow);
            b = calculateModelBMR(m, 'norms', {'L2'});
            bmr.(s{1}).kcal(i) = b.L2.kcal;
            save(ck('bmr'), 'bmr');
            fprintf('Part 7 - BMR %s %02d: %.0f kcal (control %.0f)\n', s{1}, i, b.L2.kcal, bmr.(s{1}).ctrl);
        end
    end
    save(ck('bmr'), 'bmr');
else
    load(ck('bmr'), 'bmr');
end
if stopAfter < 8, return, end

%% ==========================================================================
%  PART 8. FIGURES + TABLES
%  ==========================================================================
col = struct('male', [0.13 0.40 0.67], 'female', [0.80 0.25 0.33]);
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

% --- Fig D: read-out FVA range relative to control ----------------------------------
if exist('fva', 'var')
nR = size(readout, 1);
figure('Position', [100 100 900 500]);
W = [];  G = [];  Sx = [];
for k = 1:2
    s = sexes{k};
    wCtrl = fva.(s).ctrlMax - fva.(s).ctrlMin;
    wInd  = (fva.(s).max - fva.(s).min) ./ wCtrl;
    keep = all(~isnan(wInd), 2);
    W = [W; reshape(wInd(keep, :), [], 1)]; %#ok<AGROW>
    G = [G; reshape(repmat(1:nR, sum(keep), 1), [], 1)]; %#ok<AGROW>
    Sx = [Sx; repmat(k, sum(keep) * nR, 1)]; %#ok<AGROW>
    say('D %s: %d models; median range relative to control per read-out: %s\n', s, sum(keep), ...
        strjoin(compose('%s %.2f', string(readout(:, 2)), median(wInd(keep, :), 1)'), '; '));
end
boxchart(G, W, 'GroupByColor', categorical(Sx, [1 2], {'Male','Female'}));
set(gca, 'XTick', 1:nR, 'XTickLabel', readout(:, 2), 'XTickLabelRotation', 40);
yline(1, 'k--'); ylabel('Feasible flux range / control range');
legend('Location', 'northoutside', 'Orientation', 'horizontal');
title(sprintf('Read-out flux ranges after personalisation (window \\pm%d%%)', round(100*refWindow)));
saveFig(gcf, figPath, '', 'metabD_fvaRange');                      % not in the chapter

% --- Fig E: which input metabolite drives which read-out (Spearman, BH-FDR) ---------
% read-out value = the bound in the physiological direction of the benchmark:
% maximal release (FVA max) for release reactions, maximal uptake (FVA min) otherwise
isRelease = ismember(readout(:, 1), {'Liver_EX_glc_D(e)_[bc]','Kidney_EX_glc_D(e)_[bc]', ...
    'Adipocytes_EX_glyc(e)_[bc]','Liver_EX_bhb(e)_[bc]','Muscle_EX_ala_L(e)_[bc]','EX_urea[u]'})';
figure('Position', [100 100 1100 650]);
tl = tiledlayout(1, 2, 'TileSpacing', 'compact');
corrOut = [];
for k = 1:2
    s = sexes{k};
    Y = fva.(s).min;  Y(:, isRelease) = fva.(s).max(:, isRelease);
    keep = all(~isnan(Y), 2);
    X = conc.(s)(keep, :);  Y = Y(keep, :);
    [rho, p] = spearman(X, Y);
    rho(:, std(Y, 0, 1) < 1e-9) = NaN;  p(:, std(Y, 0, 1) < 1e-9) = NaN;   % constant read-out
    q = bhFDR(p);
    nexttile;
    imagesc(rho, 'AlphaData', ~isnan(rho)); caxis([-1 1]); colormap(divergingMap());
    set(gca, 'YTick', 1:nMet, 'YTickLabel', panel.name, 'XTick', 1:nR, 'XTickLabel', readout(:, 2), ...
        'XTickLabelRotation', 45, 'FontSize', 8);
    [ri, ci] = find(q < 0.05);
    hold on; text(ci, ri, '*', 'HorizontalAlignment', 'center', 'FontSize', 12); hold off;
    title(sprintf('%s (n = %d)', s, sum(keep)));
    [ri, ci] = find(~isnan(q));
    Tk = table(repmat(string(s), numel(ri), 1), string(metIDs(ri))', ...
        string(readout(ci, 2)), rho(sub2ind(size(rho), ri, ci)), p(sub2ind(size(p), ri, ci)), ...
        q(sub2ind(size(q), ri, ci)), 'VariableNames', {'sex','met','readout','rho','p','q'});
    if isempty(corrOut), corrOut = Tk; else, corrOut = [corrOut; Tk]; end %#ok<AGROW>
end
cb = colorbar; cb.Layout.Tile = 'east'; cb.Label.String = 'Spearman \rho';
title(tl, 'Input metabolite concentration vs read-out flux capacity (* FDR < 0.05)');
saveFig(gcf, figPath, '', 'metabE_drivers');                       % not in the chapter
writetable(sortrows(corrOut, 'q'), fullfile(resPath, 'drivers_spearman.csv'));
sig = corrOut(corrOut.q < 0.05, :);
say('E %d metabolite x read-out associations at FDR < 0.05\n', height(sig));
for r = 1:min(height(sig), 25)
    say('   %s %-8s %-24s rho %+.2f q %.2g\n', sig.sex(r), sig.met(r), sig.readout(r), sig.rho(r), sig.q(r));
end

end   % Fig D/E need Part 6

% --- Fig F: BMR ---------------------------------------------------------------------
if exist('bmr', 'var')
figure('Position', [100 100 500 420]);
hold on
for k = 1:2
    s = sexes{k};  y = bmr.(s).kcal(~isnan(bmr.(s).kcal));
    swarmX = k + 0.25 * (rand(size(y)) - 0.5);
    scatter(swarmX, y, 25, col.(s), 'filled', 'MarkerFaceAlpha', 0.6);
    plot(k + [-0.3 0.3], bmr.(s).ctrl * [1 1], 'k-', 'LineWidth', 2);
    say('F %s: control L2 BMR %.0f kcal; personalised median %.0f (range %.0f-%.0f, n = %d)\n', ...
        s, bmr.(s).ctrl, median(y), min(y), max(y), numel(y));
end
hold off
set(gca, 'XTick', 1:2, 'XTickLabel', {'Male','Female'}, 'XLim', [0.4 2.6]);
ylabel('BMR (kcal/day, L2)'); title('BMR of metabolomics-personalised WBMs (bar = control)');
saveFig(gcf, figPath, '', 'metabF_bmr');                           % not in the chapter
end   % Fig F needs Part 7

fclose(sumFile);
fprintf('Part 8 - figures in %s and %s\n', figPath, thesisFig);


%% ==========================================================================
%  LOCAL FUNCTIONS
%  ==========================================================================
function s = fmtNum(x)
% Number for a LaTeX table: 0 decimals from 100, 1 from 10, else 2.
if ~isfinite(x), s = '--';
elseif abs(x) >= 100, s = sprintf('%.0f', x);
elseif abs(x) >= 10,  s = sprintf('%.1f', x);
else,                 s = sprintf('%.2f', x);
end
end

function ok = isFeasibleWBM(m, withFloors)
% Unit whole-body maintenance flux achievable (Persephone criterion), optionally with
% the resting brain/heart ATP floors of calculateModelBMR.
if withFloors
    m.lb(contains(m.rxns, 'Brain_DM_atp')) = 3500;
    m.lb(contains(m.rxns, 'Heart_DM_atp')) = 1000;
end
ob = strcmp(m.rxns, 'Whole_body_objective_rxn');
m.c(:) = 0;  m.c(ob) = 1;  m.osenseStr = 'max';
sol = optimizeWBModel(m);
ok = sol.stat == 1 && sol.f >= 1 - 1e-6;
end

function [rho, p] = spearman(X, Y)
% Spearman correlation of every column of X with every column of Y (average ranks for
% ties); two-sided p from the t approximation, via betainc (no toolbox needed).
n = size(X, 1);
rX = rankCols(X);  rY = rankCols(Y);
rho = corr_(rX, rY);
t2 = rho.^2 .* (n - 2) ./ max(1 - rho.^2, eps);
p = betainc((n - 2) ./ (n - 2 + t2), (n - 2) / 2, 0.5);
end

function R = rankCols(X)
R = zeros(size(X));
for j = 1:size(X, 2)
    [s, ix] = sort(X(:, j));
    r = zeros(size(s));  i = 1;
    while i <= numel(s)
        k = i;
        while k < numel(s) && s(k + 1) == s(i), k = k + 1; end
        r(i:k) = (i + k) / 2;  i = k + 1;
    end
    R(ix, j) = r;
end
end

function C = corr_(A, B)
A = (A - mean(A, 1)) ./ std(A, 0, 1);  B = (B - mean(B, 1)) ./ std(B, 0, 1);
C = (A' * B) / (size(A, 1) - 1);
end

function q = bhFDR(p)
% Benjamini-Hochberg over all non-NaN entries.
q = nan(size(p));  v = find(~isnan(p));
[ps, o] = sort(p(v));  m = numel(ps);
qs = min(1, cummin(ps .* m ./ (1:m)', 'reverse'));
q(v(o)) = qs;
end

function cm = divergingMap()
n = 128;
blue = [0.13 0.40 0.67];  red = [0.80 0.25 0.33];
cm = [blue + (1 - blue) .* linspace(0, 1, n)'; flipud(red + (1 - red) .* linspace(0, 1, n)')];
end

function saveFig(h, figPath, thesisFig, name)
exportgraphics(h, fullfile(figPath, [name '.png']), 'Resolution', 300);
exportgraphics(h, fullfile(figPath, [name '.pdf']), 'ContentType', 'vector');
if exist(thesisFig, 'dir'), copyfile(fullfile(figPath, [name '.png']), thesisFig); end
end
