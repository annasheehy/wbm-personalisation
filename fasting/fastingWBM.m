%% fastingWBM.m
% Builds feasible fasting WBMs for Harvey (male) and Harvetta (female).
% Reference: Thiele et al. 2020, Mol Syst Biol e8982.
% Units: mmol/person/day.  Solver: Gurobi via COBRA Toolbox.
%
% PARTS
%   0  Setup (paths, solver)
%   1  Load models + physiological constraints, lift coupling constraints (checkpointed)
%   2  Fed models: the control models of the chapter
%   3  Fasted models (setFastingConstraints: stores released, diet closed except water
%      and ions, endogenous phosphate, hepatic glycogen, intestinal choline)
%   4  Fed and fasted BMR (L2, calculateModelBMR); resting brain/heart ATP limits
%   5  Fed -> fasting organ exchange fluxes; saves ckpt_fluxL2.mat (used by
%      fastingTables and fastingUrineCreatinine)
%
% The tests that set the fasted state (amino-acid release, phosphate) are in
% fastingSetupTests.m.

%% ==========================================================================
%  PART 0.  SETUP
%  ==========================================================================
P = persWBMpaths();
resPath        = fullfile(P.results, 'fasting');
rawCheckpoint  = fullfile(resPath, 'ckpt_rawModels.mat');
hmdbCheckpoint = fullfile(resPath, 'ckpt_HMDBconstrained.mat');
liftCheckpoint = fullfile(resPath, 'ckpt_liftedBase.mat');
feasCheckpoint = fullfile(resPath, 'ckpt_feasibleModels.mat');
fluxCheckpoint = fullfile(resPath, 'ckpt_fluxL2.mat');
if ~exist(resPath, 'dir'), mkdir(resPath); end

global CBT_LP_SOLVER CBT_QP_SOLVER
if isempty(CBT_LP_SOLVER)
    initCobraToolbox;          % slow — only on first run of the session
end
% Always assert the solver. The isempty guard above must gate ONLY the slow
% initCobraToolbox, NOT solver selection: otherwise a solver set by a previous
% session (e.g. ibm_cplex) silently persists and the BMR min-norm QP goes
% through CPLEX, which returns stat=3 (unconverged) with inflated fluxes.
okLP = changeCobraSolver('gurobi', 'LP', 0);
okQP = changeCobraSolver('gurobi', 'QP', 0);
if ~okLP || ~okQP
    warning('Gurobi not set (LP ok=%d, QP ok=%d) — still on %s. Check Gurobi install/licence.', ...
        okLP, okQP, CBT_LP_SOLVER);
end
fprintf('Solver: %s (LP)  %s (QP)  |  Output: %s\n', CBT_LP_SOLVER, CBT_QP_SOLVER, resPath);


%% ==========================================================================
%  PART 1.  LOAD MODELS + PHYSIOLOGICAL CONSTRAINTS  (checkpointed)
%  ==========================================================================

%% 1.1  Raw models -----------------------------------------------------------
if exist(rawCheckpoint, 'file')
    load(rawCheckpoint, 'harvey', 'harvetta');
    fprintf('1.1 — raw models from checkpoint.\n');
else
    harvey   = loadPSCMfile('Harvey');
    harvetta = loadPSCMfile('Harvetta');
    harvey.osenseStr = 'max';  harvetta.osenseStr = 'max';
    % Close liver biomass (regenerative growth, not fasting maintenance)
    harvey.lb(strcmp(harvey.rxns,     'Liver_biomass_reaction')) = 0;
    harvey.ub(strcmp(harvey.rxns,     'Liver_biomass_reaction')) = 0;
    harvetta.lb(strcmp(harvetta.rxns, 'Liver_biomass_reaction')) = 0;
    harvetta.ub(strcmp(harvetta.rxns, 'Liver_biomass_reaction')) = 0;
    save(rawCheckpoint, 'harvey', 'harvetta');
    fprintf('1.1 — raw models loaded and cached.\n');
end

%% 1.2  HMDB physiological constraints ---------------------------------------
% Organ-specific bounds from HMDB plasma concentrations x organ blood flow.
% Clears the objective vector -> reassign immediately after (Thiele 2020).
if exist(hmdbCheckpoint, 'file')
    load(hmdbCheckpoint, 'maleL', 'femaleL', 'ipMale', 'ipFemale');
    fprintf('1.2 — HMDB-constrained models from checkpoint.\n');
else
    sex = 'male';   standardPhysiolDefaultParameters; ipMale   = IndividualParameters;
    sex = 'female'; standardPhysiolDefaultParameters; ipFemale = IndividualParameters;
    maleL   = physiologicalConstraintsHMDBbased(harvey,   ipMale);
    femaleL = physiologicalConstraintsHMDBbased(harvetta, ipFemale);
    maleL   = changeObjective(maleL,   'Whole_body_objective_rxn');
    femaleL = changeObjective(femaleL, 'Whole_body_objective_rxn');
    maleL.osenseStr = 'max';  femaleL.osenseStr = 'max';
    save(hmdbCheckpoint, 'maleL', 'femaleL', 'ipMale', 'ipFemale');
    fprintf('1.2 — HMDB constraints applied and cached.\n');
end

%% 1.3  Lift the base models ONCE (cached) — the L2 QP fix, done up front -----
% liftCouplingConstraints reformulates ONLY the coupling constraints, which are
% invariant to the diet / storage / AA bound changes downstream. So lifting the
% base once means every model built from it (fed, fasting, ion-capped) is already
% lifted, and calculateModelBMR never re-lifts (its isModelLifted guard). This
% turns ~8 lifts per run into 0 (cached) — the ~2 h runtime was almost all lifting.
% NOTE: the cache depends only on the raw model's coupling, so it survives bound
% changes; delete ckpt_liftedBase.mat only if you change the RAW models (Part 1.1).
if exist(liftCheckpoint, 'file')
    load(liftCheckpoint, 'maleL', 'femaleL');
    fprintf('1.3 — lifted base models from checkpoint.\n');
else
    if ~isModelLifted(maleL),   maleL   = liftCouplingConstraints(maleL,   1e4, 0); end
    if ~isModelLifted(femaleL), femaleL = liftCouplingConstraints(femaleL, 1e4, 0); end
    save(liftCheckpoint, 'maleL', 'femaleL');
    fprintf('1.3 — base models lifted and cached.\n');
end


%% ==========================================================================
%  PART 2.  FED MODELS (the control models of the chapter)
%  ==========================================================================
EUAverageDietNew
maleFed   = setDietConstraints(setFeedingFastingConstraints(maleL,   'feeding', -1e6, 1e6), Diet);
femaleFed = setDietConstraints(setFeedingFastingConstraints(femaleL, 'feeding', -1e6, 1e6), Diet);
maleFed.osenseStr = 'max';  femaleFed.osenseStr = 'max';
sol = optimizeWBModel(maleFed);   fprintf('2 — fed male:   stat=%d  obj=%.4f\n', sol.stat, sol.f);
sol = optimizeWBModel(femaleFed); fprintf('2 — fed female: stat=%d  obj=%.4f\n', sol.stat, sol.f);


%% ==========================================================================
%  PART 3.  FASTING MODEL
%  The complete, justified fasting constraint set is encapsulated in
%  setFastingConstraints: storage release (organ sinks drain), diet closed,
%  water + dissolved ions available, endogenous phosphate (bone Pi — the model
%  has no bone organ), and sex-specific hepatic glycogen release + intestinal
%  choline from blood. The defaults reproduce the values derived across the
%  investigations; see setFastingConstraints for the per-knob rationale and the
%  (still provisional) TODOs.
%  ==========================================================================
maleFast   = setFastingConstraints(maleL,   'male');
femaleFast = setFastingConstraints(femaleL, 'female');
maleFast.osenseStr = 'max';  femaleFast.osenseStr = 'max';
sol = optimizeWBModel(maleFast);   fprintf('3 — fasting male:   stat=%d  obj=%.4f\n', sol.stat, sol.f);
sol = optimizeWBModel(femaleFast); fprintf('3 — fasting female: stat=%d  obj=%.4f\n', sol.stat, sol.f);

save(feasCheckpoint, 'harvey', 'harvetta', 'maleL', 'femaleL', ...
     'maleFed', 'femaleFed', 'maleFast', 'femaleFast');
fprintf('3 — feasible fasting models saved to checkpoint.\n');


%% ==========================================================================
%  PART 4.  BMR (L2) AND RESTING ATP LIMITS
%  BMR by the basal metabolic flux method of Thiele 2020 (calculateModelBMR): resting
%  brain/heart ATP floors, whole-body objective fixed to 1, min ||v||_2 (QP; lifted
%  coupling, barrier, crossover off), ATP consumption summed and converted with
%  64 kJ/mol. The fed and fasted flux vectors are saved for Part 5.
%  ==========================================================================
if ~exist('maleFast', 'var'), load(feasCheckpoint); end

%% 4.1  Fed and fasted BMR (L2) ----------------------------------------------
mods = {maleFed, femaleFed; maleFast, femaleFast};    % rows fed / fasted, columns male / female
state = {'fed', 'fasted'};  sexName = {'male', 'female'};
v = cell(2, 2);
for st = 1:2
    for sx = 1:2
        b = calculateModelBMR(mods{st, sx});
        v{st, sx} = b.v;
        fprintf('4.1 - %s %s BMR (L2): %.0f kcal/day (stat %d)\n', state{st}, sexName{sx}, b.kcal, b.stat);
    end
end
[fluxFedM, fluxFedF, fluxFastM, fluxFastF] = deal(v{1, 1}, v{1, 2}, v{2, 1}, v{2, 2});

%% 4.2  Resting brain/heart ATP limits of the fasted models ------------------
% The BMR solve requires Brain_DM_atp >= 3500 and Heart_DM_atp >= 1000. With the
% whole-body objective fixed at 1, the largest ATP each organ can be supplied (LP)
% shows whether a fasted model can meet these demands.
sexCeil = {'male', maleFast; 'female', femaleFast};
for c = 1:2
    m  = setBMRConditions(sexCeil{c,2}, [], []);              % objective fixed to 1, no ATP floors
    mB = m;  mB.c(:) = 0;  mB.c(contains(m.rxns, 'Brain_DM_atp')) = 1;  mB.osenseStr = 'max';
    mH = m;  mH.c(:) = 0;  mH.c(contains(m.rxns, 'Heart_DM_atp')) = 1;  mH.osenseStr = 'max';
    sB = optimizeWBModel(mB);  sH = optimizeWBModel(mH);
    fprintf('4.2 - %-6s: brain ATP max %6.0f (demand 3500), heart ATP max %6.0f (demand 1000)\n', ...
        sexCeil{c,1}, sB.f, sH.f);
end


%% ==========================================================================
%  PART 5.  FED -> FASTING ORGAN EXCHANGE FLUXES  (L2 min-norm QP)
%  Data-led survey of the inter-organ fuel economy: every organ<->blood exchange
%  reaction (<Organ>_EX_<met>(e)_[bc], and brain _[csf]) with its fed and fasting
%  L2 flux and the change. Exchange-flux SIGN is meaningful (uptake vs release), so
%  these are printed RAW for manual interpretation — no summing, no cherry-picking.
%  Console shows the global top movers per sex (small fuels, e.g. ketones, carry
%  tiny flux and won't rank there but ARE in the Excel). Excel: one sheet per organ,
%  sorted by |change| (fedFastingExchange_<sex>.xlsx) — browse / Ctrl-F glc_D, bhb...
%  Re-runnable on its own: loads the flux vectors + models from the checkpoints.
%  ==========================================================================
if ~exist('fluxFastM', 'var') && exist(fluxCheckpoint, 'file')
    load(fluxCheckpoint, 'fluxFedM', 'fluxFedF', 'fluxFastM', 'fluxFastF');
end
if ~exist('maleFast', 'var'), load(feasCheckpoint, 'maleFast', 'femaleFast'); end

if exist('fluxFedM', 'var') && exist('fluxFastM', 'var') && ~isempty(fluxFedM) && ~isempty(fluxFastM)
    save(fluxCheckpoint, 'fluxFedM', 'fluxFedF', 'fluxFastM', 'fluxFastF');   % for standalone re-runs
    topShow = 20;
    vmhAbbr = {};  vmhName = {};                     % VMH abbreviation -> full metabolite name
    try
        vdb = loadVMHDatabase;  vmhAbbr = vdb.metabolites(:,1);  vmhName = vdb.metabolites(:,2);
    catch
        warning('Part 5: loadVMHDatabase failed — using metabolite abbreviations.');
    end
    fluxSex = {'male', maleFast, fluxFedM, fluxFastM; 'female', femaleFast, fluxFedF, fluxFastF};
    for q = 1:2
        m = fluxSex{q,2};  vFed = fluxSex{q,3};  vFast = fluxSex{q,4};
        if isempty(vFed) || isempty(vFast)
            fprintf('\n5 — %s: L2 flux missing (fed or fasting did not solve) — skipped.\n', fluxSex{q,1});
            continue
        end
        % Organ<->blood exchange reactions: <Organ>_EX_<met>(e)_[bc] (brain _[csf]).
        isExch = ~cellfun(@isempty, regexp(m.rxns, '^[A-Za-z0-9]+_EX_.*_\[(bc|csf)\]$', 'once'));
        idx    = find(isExch);
        orgTok = regexp(m.rxns(idx), '^([A-Za-z0-9]+)_EX_', 'tokens', 'once');
        organ  = cellfun(@(c) c{1}, orgTok, 'UniformOutput', false);
        fed    = full(vFed(idx));  fast = full(vFast(idx));  d = fast - fed;

        % Per exchange: the blood-side metabolite, its VMH full name, and the
        % uptake/release direction DERIVED FROM STOICHIOMETRY (sign of flux x the
        % blood-metabolite coefficient) — correct regardless of how the reaction is
        % written. release = organ -> blood; uptake = blood -> organ.
        nE = numel(idx);  metName = repmat({''}, nE, 1);
        fedDir = repmat({'-'}, nE, 1);  fastDir = repmat({'-'}, nE, 1);
        for e = 1:nE
            mr = find(m.S(:, idx(e)) ~= 0);
            bc = mr(find(~cellfun(@isempty, regexp(m.mets(mr), '\[(bc|csf)\]$', 'once')), 1));
            if isempty(bc), continue; end
            abbr = regexprep(m.mets{bc}, '\[(bc|csf)\]$', '');
            hit  = find(strcmp(vmhAbbr, abbr), 1);
            if ~isempty(hit), metName{e} = vmhName{hit}; else, metName{e} = abbr; end
            coeff = m.S(bc, idx(e));
            if     fed(e)*coeff  >  1e-9, fedDir{e}  = 'release'; elseif fed(e)*coeff  < -1e-9, fedDir{e}  = 'uptake'; end
            if     fast(e)*coeff >  1e-9, fastDir{e} = 'release'; elseif fast(e)*coeff < -1e-9, fastDir{e} = 'uptake'; end
        end

        % Excel: one sheet per organ, sorted by |change|.
        xlsFile = fullfile(resPath, sprintf('fedFastingExchange_%s.xlsx', fluxSex{q,1}));
        if exist(xlsFile, 'file'), delete(xlsFile); end
        organs = unique(organ);
        for oi = 1:numel(organs)
            sel = strcmp(organ, organs{oi});
            T = table(m.rxns(idx(sel)), metName(sel), fed(sel), fast(sel), d(sel), fedDir(sel), fastDir(sel), abs(d(sel)), ...
                'VariableNames', {'reaction','metabolite','fed','fasting','delta','fedDir','fastDir','absDelta'});
            T = sortrows(T, 'absDelta', 'descend');  T.absDelta = [];
            writetable(T, xlsFile, 'Sheet', organs{oi});
        end

        % Console: top movers, water/O2/CO2/HCO3 hidden so fuels surface (Excel has all).
        keep  = cellfun(@isempty, regexp(m.rxns(idx), '_(h2o|o2|co2|hco3)[\(\[]', 'once')) & abs(d) > 1e-3;
        cShow = find(keep);
        [~, o] = sort(abs(d(cShow)), 'descend');
        fprintf('\n7 — %s: organ<->blood exchange, top %d fed->fasting changes (water/O2/CO2 hidden; all in %s):\n', ...
            fluxSex{q,1}, topShow, xlsFile);
        fprintf('   %-28s %-20s %9s %9s %9s  %s\n', 'reaction', 'metabolite', 'fed', 'fasting', 'delta', 'fastDir');
        for k = 1:min(topShow, numel(o))
            j = cShow(o(k));  nm = metName{j};
            fprintf('   %-28s %-20s %9.1f %9.1f %9.1f  %s\n', ...
                m.rxns{idx(j)}, nm(1:min(numel(nm),20)), fed(j), fast(j), d(j), fastDir{j});
        end
    end
else
    warning('Part 5: L2 flux vectors not available — run Part 4 first.');
end
