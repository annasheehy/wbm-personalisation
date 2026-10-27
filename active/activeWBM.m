%% activeWBM.m
% Thesis Ch5, Results 2.4 (minimal): personalising physiological state to physical activity.
%
% Control: the fed Harvey/Harvetta of fastingWBM.m (diet is the only energy source; in the
% feeding state storage sinks can fill but not release).
%   1. Activity capacity: maximal extra muscle ATP hydrolysis (Muscle_DM_atp) under BMR
%      conditions (whole-body maintenance = 1, brain/heart ATP floors), in kcal/day, with
%      the reference O2 uptake bound and with it doubled (does O2 or diet energy limit?).
%      Compared with the ~500 kcal/day left for activity in Thiele et al. 2020.
%   2. One running bout: extra muscle ATP = (METs - 1) x body weight x hours (1 MET =
%      1 kcal/kg/h), imposed as a floor on Muscle_DM_atp; L2 (min-norm) flux as in
%      calculateModelBMR. Scenario 'diet' as the control; scenario 'stores' additionally
%      lets liver/muscle glycogen and adipose triacylglycerol be released; 'o2' raises the
%      whole-body O2 uptake bound by the extra O2 of the bout (from the METs); 'o2muscle'
%      also raises the muscle O2 uptake bound by the same amount (the resting muscle O2
%      bound of Harvey/Harvetta caps capacity otherwise), i.e. the extra O2 is used by the
%      working muscle. Each scenario's capacity (LP) is recorded; the bout is solved only
%      where it fits.
%      Read-outs: model energy expenditure vs BMR + activity cost, RQ, organ fuel exchanges.
%
% Energy: kcal = ATP (mmol/day) / 1000 * 64 kJ/mol * 0.239006 kcal/kJ (as calculateModelBMR).
% Anna Sheehy 2026.

%% Setup ------------------------------------------------------------------------------
runMETs    = 9.3;     % 2024 Adult Compendium (Herrmann et al. 2024), code 12050: running 6-6.3 mph (10 min/mile)
runMinutes = 60;
kcalPerMol = 64 * 0.239006;
storeRxns  = {'Liver_sink_glygn2(c)', 'Muscle_sink_glygn2(c)', 'Adipocytes_sink_tag_hs(c)'};
fuel = {'EX_o2[a]', 'O2 uptake'; 'EX_co2[a]', 'CO2 output'; ...
        'Muscle_EX_glc_D(e)_[bc]', 'Muscle glucose'; 'Muscle_EX_lac_L(e)_[bc]', 'Muscle lactate'; ...
        'Liver_EX_glc_D(e)_[bc]', 'Liver glucose'; 'Liver_EX_lac_L(e)_[bc]', 'Liver lactate'; ...
        'Adipocytes_EX_glyc(e)_[bc]', 'Adipose glycerol'; 'Muscle_EX_ala_L(e)_[bc]', 'Muscle alanine'};

P = persWBMpaths();
fastPath  = fullfile(P.results, 'fasting');
resPath   = fullfile(P.results, 'active');
if ~exist(resPath, 'dir'), mkdir(resPath); end

if exist('startDelay', 'var'), pause(startDelay); end
global CBT_LP_SOLVER
if isempty(CBT_LP_SOLVER), initCobraToolbox(false); end
changeCobraSolver('gurobi', 'LP', 0);
changeCobraSolver('gurobi', 'QP', 0);

S = load(fullfile(fastPath, 'ckpt_feasibleModels.mat'), 'maleFed', 'femaleFed');
IPd = load(fullfile(fastPath, 'ckpt_HMDBconstrained.mat'), 'ipMale', 'ipFemale');
ctrl = struct('male', S.maleFed, 'female', S.femaleFed);
ipar = struct('male', IPd.ipMale, 'female', IPd.ipFemale);
clear S IPd
sexes = {'male', 'female'};

res = struct('male', struct(), 'female', struct());
if exist(fullfile(resPath, 'activeResults.mat'), 'file')     % resume: keep finished scenarios
    R0 = load(fullfile(resPath, 'activeResults.mat'), 'res');  res = R0.res;  clear R0
end
for k = 1:2
    s = sexes{k};  m0 = ctrl.(s);
    mDM = find(contains(m0.rxns, 'Muscle_DM_atp'));
    if numel(mDM) ~= 1
        error('Expected one Muscle_DM_atp reaction, found %d: %s', numel(mDM), strjoin(m0.rxns(mDM), ', '));
    end
    o2 = find(strcmp(m0.rxns, 'EX_o2[a]'));
    mO2 = find(strcmp(m0.rxns, 'Muscle_EX_o2(e)_[bc]'));
    missing = setdiff([fuel(:, 1); storeRxns(:)], m0.rxns);
    if ~isempty(missing), warning('activeWBM:missing', 'Not in %s model: %s', s, strjoin(missing, ', ')); end

    % 1. capacity under BMR conditions (as calculateModelBMR / fastingTargetFVA)
    cap = nan(1, 2);
    for o2scale = [1 2]
        m = m0;  m.lb(o2) = o2scale * m0.lb(o2);
        cap(o2scale) = maxMuscleATP(m, mDM);
    end
    res.(s).capacityATP  = cap;                       % mmol/day [reference O2, O2 x2]
    res.(s).capacityKcal = cap / 1000 * kcalPerMol;
    BW = ipar.(s).bodyWeight;
    res.(s).capacityRunMin = res.(s).capacityKcal / ((runMETs - 1) * BW) * 60;
    fprintf('%s: activity capacity %.0f kcal/day (O2 x2: %.0f) = %.0f min running at %.1f METs\n', s, ...
        res.(s).capacityKcal(1), res.(s).capacityKcal(2), res.(s).capacityRunMin(1), runMETs);

    % 2. running bout
    boutKcal = (runMETs - 1) * BW * runMinutes / 60;
    boutATP  = boutKcal / kcalPerMol * 1000;
    res.(s).boutKcal = boutKcal;  res.(s).boutATP = boutATP;
    % extra O2 for the bout: 1 MET = 3.5 mL O2/kg/min (Compendium); 22.4 mL per mmol
    extraO2 = (runMETs - 1) * 3.5 * BW * runMinutes / 22.4;
    res.(s).extraO2 = extraO2;
    if ~(isfield(res.(s), 'rest') && ~isempty(res.(s).rest.fluxes))
        rest = calculateModelBMR(m0, 'norms', {'L2'});
        res.(s).rest = struct('kcal', rest.L2.kcal, 'fluxes', pickFluxes(rest.L2.v, m0, fuel));
    end
    for scen = {'diet', 'stores', 'o2', 'o2stores', 'o2muscle', 'o2musclestores'}
        m = m0;
        if contains(scen{1}, 'stores')
            st = ismember(m.rxns, storeRxns);  m.lb(st) = -1e6;
        end
        if startsWith(scen{1}, 'o2')
            m.lb(o2) = m0.lb(o2) - extraO2;               % uptake is negative
        end
        if contains(scen{1}, 'muscle')
            m.lb(mO2) = m0.lb(mO2) - extraO2;             % extra O2 taken up by the working muscle
        end
        capKcal = maxMuscleATP(m, mDM) / 1000 * kcalPerMol;
        if isfield(res.(s), scen{1}) && ~isempty(res.(s).(scen{1}).fluxes)   % bout solved earlier
            res.(s).(scen{1}).capacityKcal = capKcal;
            continue
        end
        if ~(capKcal >= boutKcal)                         % bout does not fit (or capacity LP failed)
            res.(s).(scen{1}) = struct('stat', 0, 'kcal', NaN, 'fluxes', [], 'capacityKcal', capKcal);
            fprintf('%s bout (%s): capacity %.0f kcal < bout %.0f kcal, not solved\n', s, scen{1}, capKcal, boutKcal);
            continue
        end
        m.lb(mDM) = boutATP;
        a = calculateModelBMR(m, 'norms', {'L2'});
        r = struct('stat', a.L2.stat, 'kcal', a.L2.kcal, 'fluxes', pickFluxes(a.L2.v, m, fuel), ...
                   'capacityKcal', capKcal);
        if contains(scen{1}, 'stores') && ~isempty(a.L2.v)
            r.storeRelease = a.L2.v(ismember(m.rxns, storeRxns))';
        end
        res.(s).(scen{1}) = r;
        fprintf('%s bout (%s): stat %d, model EE %.0f kcal vs rest %.0f + bout %.0f = %.0f\n', s, scen{1}, ...
            r.stat, r.kcal, res.(s).rest.kcal, boutKcal, res.(s).rest.kcal + boutKcal);
        save(fullfile(resPath, 'activeResults.mat'), 'res', 'runMETs', 'runMinutes', 'fuel', 'storeRxns');
    end
    save(fullfile(resPath, 'activeResults.mat'), 'res', 'runMETs', 'runMinutes', 'fuel', 'storeRxns');
end

% summary table
rows = {};
for k = 1:2
    s = sexes{k};
    for scen = {'rest', 'diet', 'stores', 'o2', 'o2stores', 'o2muscle', 'o2musclestores'}
        if ~isfield(res.(s), scen{1}), continue, end
        r = res.(s).(scen{1});
        cap = NaN;  if isfield(r, 'capacityKcal'), cap = r.capacityKcal; end
        if isempty(r.fluxes)
            rows(end+1, :) = [{s, scen{1}, cap, NaN, NaN}, num2cell(nan(1, size(fuel, 1)))]; %#ok<SAGROW>
            continue
        end
        rq = r.fluxes(2) / -r.fluxes(1);
        rows(end+1, :) = [{s, scen{1}, cap, r.kcal, rq}, num2cell(r.fluxes)]; %#ok<SAGROW>
    end
end
T = cell2table(rows, 'VariableNames', [{'sex', 'state', 'capacityKcal', 'EEkcal', 'RQ'}, matlab.lang.makeValidName(fuel(:, 2))']);
writetable(T, fullfile(resPath, 'activeSummary.csv'));
disp(T)

% LaTeX table: rest, run cost, and capacity, energy expenditure and RQ per scenario and sex
activeStateTable(res, runMETs, fullfile(P.tables, 'activeState.tex'));

function f = maxMuscleATP(m, mDM)
% maximal extra muscle ATP hydrolysis (mmol/day) under BMR conditions; NaN if not solved
m.lb(contains(m.rxns, 'Brain_DM_atp')) = 3500;
m.lb(contains(m.rxns, 'Heart_DM_atp')) = 1000;
ob = strcmp(m.rxns, 'Whole_body_objective_rxn');  m.lb(ob) = 1;  m.ub(ob) = 1;
m.c(:) = 0;  m.c(mDM) = 1;  m.osenseStr = 'max';
sol = optimizeWBModel(m);
f = NaN;  if sol.stat == 1, f = sol.f; end
end

function f = pickFluxes(v, m, fuel)
f = [];
if isempty(v), return, end
[ok, idx] = ismember(fuel(:, 1), m.rxns);
f = nan(1, size(fuel, 1));  f(ok) = v(idx(ok))';
end
