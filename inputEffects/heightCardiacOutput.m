%% heightCardiacOutput.m
% Thesis Ch5, Results 2.1: height through cardiac output and blood volume, and haematocrit
% through red cell mass.
% In the default parameters, cardiac output (CO) = heart rate x stroke volume, and height
% changes only the blood volume, which no constraint uses; calcOrganFract (PSCM toolbox) also
% replaces the personalised blood volume by an estimate from body weight when it sets the
% blood cell weights of the whole-body maintenance reaction. This script tests:
%   - CO from height: individualizedLabReport optionCardiacOutput = 4, CO from the Du Bois body
%     surface area and 125 mL O2/min per m2 (Fick principle), at low, default and high height,
%     converted to mL/min (see coOption4);
%   - the same with the blood cell weights set from the personalised blood volume
%     (calcOrganFractPersonal), with red cell weight scaled by haematocrit / default;
%   - haematocrit at its low and high healthy value with these blood cell weights.
% References: the control (oatFVA/oatL2 control jobs of metabolomicsOAT.m), and the default
% parameters with organ weights from calcOrganFract ('std') or calcOrganFractPersonal ('pers').
% Read-outs: FVA of the 11 organ exchanges and QP-FBA BMR and O2 uptake.
%
% Anna Sheehy 2026.

if ~exist('sexList', 'var'), sexList = {'male', 'female'}; end
% CO (mL/min) = 125 mL O2/min per m2 x BSA / (200 - 150 mL O2 per L blood) x 1000 mL/L: option 4
% of individualizedLabReport, with its conversion to mL/min (fixed in the COBRA Toolbox; before
% the fix it multiplied by 60*24, giving L/day, 1.44 times the mL/min the constraints use).
coOption4 = @(ip) 125 * bsaDuBois(ip.bodyWeight, ip.Height) / (200 - 150) * 1000;
readout = oatReadouts();  readout = readout(:, 1);

P = persWBMpaths();
fastPath = fullfile(P.results, 'fasting');
resPath  = fullfile(P.results, 'metabolomics');
outDir   = fullfile(resPath, 'heightCO');  if ~exist(outDir, 'dir'), mkdir(outDir); end
if exist('startDelay', 'var'), pause(startDelay); end
global CBT_LP_SOLVER
if isempty(CBT_LP_SOLVER), initCobraToolbox(false); end
changeCobraSolver('gurobi', 'LP', 0);  changeCobraSolver('gurobi', 'QP', 0);

S  = load(fullfile(fastPath, 'ckpt_feasibleModels.mat'), 'maleFed', 'femaleFed');
IPd = load(fullfile(fastPath, 'ckpt_HMDBconstrained.mat'), 'ipMale', 'ipFemale');
Rw = load(fullfile(fastPath, 'ckpt_rawModels.mat'), 'harvey', 'harvetta');
ctrl = struct('male', S.maleFed, 'female', S.femaleFed);
ipar = struct('male', IPd.ipMale, 'female', IPd.ipFemale);
raw  = struct('male', Rw.harvey, 'female', Rw.harvetta);
clear S IPd Rw
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
        m = onLiftedModel(ctrl.(sx), buildFedWBM(raw.(sx), ip, Diet, cases{c, 5}, ip0.Hematocrit));
        F = fastingTargetFVA(m, readout);
        b = calculateModelBMR(m);
        o2 = find(strcmp(m.rxns, 'EX_o2[a]'));
        R = struct('case', cases{c, 1}, 'height', cases{c, 2}, 'hct', cases{c, 3}, 'coMode', cases{c, 4}, ...
            'orgMode', cases{c, 5}, 'CO', ip.CardiacOutput, 'COdefault', ip0.CardiacOutput, 'BV', ip.BloodVolume, ...
            'F', F, 'kcal', b.kcal, 'stat', b.stat, 'o2', b.v(o2), 'v', single(b.v), 'secs', toc(t0));
        save(f, '-struct', 'R');
        fprintf('%s %s: height %.1f cm, Hct %.3f, CO %.0f mL/min, BV %.0f mL; BMR %.0f kcal (stat %d), O2 %.0f; %.0f s\n', ...
            sx, cases{c, 1}, R.height, R.hct, R.CO, R.BV, R.kcal, R.stat, R.o2, R.secs);
    end
end
