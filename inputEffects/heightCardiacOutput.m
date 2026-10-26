%% heightCardiacOutput.m
% Thesis Ch5, Results 2.1: height through cardiac output and blood volume, and haematocrit
% through red cell mass.
% In the default parameters, cardiac output (CO) = heart rate x stroke volume, and height
% changes only the blood volume, which no constraint uses; calcOrganFract (PSCM toolbox) also
% replaces the personalised blood volume by an estimate from body weight when it sets the
% blood cell weights of the whole-body maintenance reaction. This script tests:
%   - CO from height: individualizedLabReport optionCardiacOutput = 4, CO from the Du Bois body
%     surface area and 125 mL O2/min per m2 (Fick principle), at low, default and high height;
%   - the same with the blood cell weights set from the personalised blood volume
%     (calcOrganFractPersonal), with red cell weight scaled by haematocrit / default;
%   - haematocrit at its low and high healthy value with these blood cell weights.
% References: the control (oatFVA/oatL2 control jobs of metabolomicsOAT.m), and the default
% parameters with organ weights from calcOrganFract ('std') or calcOrganFractPersonal ('pers').
% Read-outs: FVA of the 11 organ exchanges and QP-FBA BMR and O2 uptake.
%
% Anna Sheehy 2026.

if ~exist('sexList', 'var'), sexList = {'male', 'female'}; end
objRxn = 'Whole_body_objective_rxn';
bsaDuBois = @(W, H) 0.007184 * W^0.425 * H^0.725;                      % m2; W kg, H cm
coOption4 = @(ip) 0.125 * 1000 * bsaDuBois(ip.bodyWeight, ip.Height) * 60 * 24 / (200 - 150);  % as individualizedLabReport
readout = {'Liver_EX_glc_D(e)_[bc]', 'Brain_EX_glc_D(e)_[csf]', 'Muscle_EX_glc_D(e)_[bc]', 'RBC_EX_glc_D(e)_[bc]', ...
           'Kidney_EX_glc_D(e)_[bc]', 'Adipocytes_EX_glyc(e)_[bc]', 'Liver_EX_bhb(e)_[bc]', 'Brain_EX_bhb(e)_[csf]', ...
           'Muscle_EX_ala_L(e)_[bc]', 'Liver_EX_ala_L(e)_[bc]', 'EX_urea[u]'};

P = persWBMpaths();
fastPath = fullfile(P.results, 'fasting');
resPath  = fullfile(P.results, 'metabolomics');
outDir   = fullfile(resPath, 'heightCO');  if ~exist(outDir, 'dir'), mkdir(outDir); end
if exist('startDelay', 'var'), pause(startDelay); end
global CBT_LP_SOLVER
if isempty(CBT_LP_SOLVER), initCobraToolbox(false); end
changeCobraSolver('gurobi', 'LP', 0);  changeCobraSolver('gurobi', 'QP', 0);

S  = load(fullfile(fastPath, 'ckpt_feasibleModels.mat'), 'maleFed', 'femaleFed');
P  = load(fullfile(fastPath, 'ckpt_HMDBconstrained.mat'), 'ipMale', 'ipFemale');
Rw = load(fullfile(fastPath, 'ckpt_rawModels.mat'), 'harvey', 'harvetta');
ctrl = struct('male', S.maleFed, 'female', S.femaleFed);
ipar = struct('male', P.ipMale, 'female', P.ipFemale);
raw  = struct('male', Rw.harvey, 'female', Rw.harvetta);
clear S P Rw
EUAverageDietNew;                                                     % defines Diet
lims = healthyPhysiologyRanges(fullfile(P.references, 'healthyPhysiologyRanges.csv'), ipar);

for s = sexList
    sx = s{1};  ip0 = ipar.(sx);  H = lims.(sx).Height;  Hct = lims.(sx).Hematocrit;
    % name, height, haematocrit, CO ('default' | 'opt4'), organ weights ('none' | 'std' | 'pers')
    cases = {'defaultHeight',     ip0.Height, ip0.Hematocrit, 'opt4',    'none'
             'lowHeight',         H(1),       ip0.Hematocrit, 'opt4',    'none'
             'highHeight',        H(2),       ip0.Hematocrit, 'opt4',    'none'
             'std_default',       ip0.Height, ip0.Hematocrit, 'default', 'std'
             'pers_default',      ip0.Height, ip0.Hematocrit, 'default', 'pers'
             'opt4pers_default',  ip0.Height, ip0.Hematocrit, 'opt4',    'pers'
             'opt4pers_lowH',     H(1),       ip0.Hematocrit, 'opt4',    'pers'
             'opt4pers_highH',    H(2),       ip0.Hematocrit, 'opt4',    'pers'
             'pers_lowHct',       ip0.Height, Hct(1),         'default', 'pers'
             'pers_highHct',      ip0.Height, Hct(2),         'default', 'pers'};
    for c = 1:size(cases, 1)
        f = fullfile(outDir, sprintf('%s_%s.mat', sx, cases{c, 1}));
        if exist(f, 'file'), continue, end
        if exist('skipCases', 'var') && ismember(cases{c, 1}, skipCases), continue, end   % run by another session
        ip = recalcDependentParams(ip0, 'Height', cases{c, 2});
        ip.Hematocrit = cases{c, 3};
        if strcmp(cases{c, 4}, 'opt4'), ip.CardiacOutput = coOption4(ip); end
        t0 = tic;
        mm = buildPhysModel(raw.(sx), ip, cases{c, 5}, objRxn, Diet, ip0.Hematocrit);
        m = onLiftedControl(ctrl.(sx), mm);
        F = fastingTargetFVA(m, readout);
        b = calculateModelBMR(m, 'norms', {'L2'});
        o2 = find(strcmp(m.rxns, 'EX_o2[a]'));
        R = struct('case', cases{c, 1}, 'height', cases{c, 2}, 'hct', cases{c, 3}, 'coMode', cases{c, 4}, ...
            'orgMode', cases{c, 5}, 'CO', ip.CardiacOutput, 'COdefault', ip0.CardiacOutput, 'BV', ip.BloodVolume, ...
            'F', F, 'kcal', b.L2.kcal, 'stat', b.L2.stat, 'o2', b.L2.v(o2), 'v', single(b.L2.v), 'secs', toc(t0));
        save(f, '-struct', 'R');
        fprintf('%s %s: height %.1f cm, Hct %.3f, CO %.0f mL/min, BV %.0f mL; BMR %.0f kcal (stat %d), O2 %.0f; %.0f s\n', ...
            sx, cases{c, 1}, R.height, R.hct, R.CO, R.BV, R.kcal, R.stat, R.o2, R.secs);
    end
end

function m = buildPhysModel(rawModel, ip, orgMode, objRxn, Diet, hctRef)
% as in metabolomicsOAT.m; orgMode 'std' recomputes organ weights with calcOrganFract,
% 'pers' with calcOrganFractPersonal (personalised blood volume, red cells scaled by Hct)
m = rawModel;
if ~strcmp(orgMode, 'none')
    if strcmp(orgMode, 'pers')
        [org, ~, owf, ip] = calcOrganFractPersonal(m, ip, hctRef);
    else
        [org, ~, owf, ip] = calcOrganFract(m, ip);
    end
    m = adjustWholeBodyRxnCoeff(m, org, owf);
end
m = physiologicalConstraintsHMDBbased(m, ip);
m = changeObjective(m, objRxn);  m.osenseStr = 'max';
m = setDietConstraints(setFeedingFastingConstraints(m, 'feeding', -1e6, 1e6), Diet);
m.osenseStr = 'max';
end

function m = onLiftedControl(liftedCtrl, mm)
% as in metabolomicsOAT.m: put the rebuilt bounds onto the lifted control
assert(isequal(liftedCtrl.rxns, mm.rxns) && isequal(liftedCtrl.mets, mm.mets), ...
    'onLiftedControl: reaction or metabolite lists differ');
m = liftedCtrl;
m.lb = mm.lb;  m.ub = mm.ub;  m.c = mm.c;
d = find(any(m.S ~= mm.S, 1));
if numel(d) > 1
    error('onLiftedControl: %d stoichiometric columns differ (expected 0 or 1)', numel(d));
end
m.S(:, d) = mm.S(:, d);
end
