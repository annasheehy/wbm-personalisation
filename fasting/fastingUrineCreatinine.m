%% fastingUrineCreatinine.m
% Thesis Ch5, fasted state: urinary urea excretion of the fasted models with a measured urine
% creatinine concentration. Urinary limits given per mmol creatinine scale with the urine
% creatinine concentration, whose default (0.5-1.2 mg/dL, standardPhysiolDefaultParameters) is
% the serum range. The fasted models of fastingWBM.m are rebuilt with urine creatinine at the
% low and high end of the range tested in metabolomicsOAT.m (references/
% healthyPhysiologyRanges.csv), and feasibility, BMR (QP-FBA) and the urinary urea flux
% (QP-FBA and FVA) are compared with the models as built.
% OUTPUT: results/fasting/urineCreatinine/<sex>_<level>.mat and summary printed
%
% Anna Sheehy 2026.

P = persWBMpaths();
resPath = fullfile(P.results, 'fasting');
outPath = fullfile(resPath, 'urineCreatinine');
if ~exist(outPath, 'dir'), mkdir(outPath); end
global CBT_LP_SOLVER
if isempty(CBT_LP_SOLVER), initCobraToolbox(false); end
changeCobraSolver('gurobi', 'LP', 0);  changeCobraSolver('gurobi', 'QP', 0);

Rw  = load(fullfile(resPath, 'ckpt_rawModels.mat'), 'harvey', 'harvetta');
Lf  = load(fullfile(resPath, 'ckpt_liftedBase.mat'), 'maleL', 'femaleL');
Fm  = load(fullfile(resPath, 'ckpt_feasibleModels.mat'), 'maleFast', 'femaleFast');
IPd = load(fullfile(resPath, 'ckpt_HMDBconstrained.mat'), 'ipMale', 'ipFemale');
L2  = load(fullfile(resPath, 'ckpt_fluxL2.mat'), 'fluxFastM', 'fluxFastF');
ipar = struct('male', IPd.ipMale, 'female', IPd.ipFemale);
lims = healthyPhysiologyRanges(fullfile(P.references, 'healthyPhysiologyRanges.csv'), ipar);
cases = {'male', Rw.harvey, Lf.maleL, Fm.maleFast, L2.fluxFastM; 'female', Rw.harvetta, Lf.femaleL, Fm.femaleFast, L2.fluxFastF};
urea = 'EX_urea[u]';

for k = 1:2
    [s, raw, lifted, asBuilt, vAsBuilt] = cases{k, :};
    i = strcmp(asBuilt.rxns, urea);
    m0 = lifted;  mm = physiologicalConstraintsHMDBbased(raw, ipar.(s));  m0.lb = mm.lb;  m0.ub = mm.ub;
    m0 = setFastingConstraints(m0, s);                            % check: default parameters rebuild the as-built model
    fprintf('%s rebuild check: max bound difference to the as-built fasted model %.3g\n', s, ...
        max([abs(m0.lb - asBuilt.lb); abs(m0.ub - asBuilt.ub)]));
    F0 = fastingTargetFVA(asBuilt, {urea});
    fprintf('\n==== %s as built: urine creatinine %.1f-%.1f mg/dL; urinary urea QP-FBA %.2f, FVA [%.2f %.2f]\n', s, ...
        ipar.(s).MConUrCreatinineMin, ipar.(s).MConUrCreatinineMax, full(vAsBuilt(i)), F0.fvaMin, F0.fvaMax);
    for L = 1:2
        f = fullfile(outPath, sprintf('%s_%d.mat', s, L));
        if exist(f, 'file'), load(f, 'r'); else
            uc = lims.(s).UrineCreatinine(L);
            mm = physiologicalConstraintsHMDBbased(raw, recalcDependentParams(ipar.(s), 'UrineCreatinine', uc));
            m = lifted;  m.lb = mm.lb;  m.ub = mm.ub;               % HMDB bounds on the lifted base, as fastingWBM.m
            m = setFastingConstraints(m, s);
            m = changeObjective(m, 'Whole_body_objective_rxn');  m.osenseStr = 'max';
            sol = optimizeWBModel(m);
            r = struct('uc', uc, 'feasStat', sol.stat, 'kcal', NaN, 'ureaQP', NaN, 'ureaFVA', [NaN NaN]);
            if sol.stat == 1
                b = calculateModelBMR(m);
                r.kcal = b.kcal;  r.ureaQP = b.v(strcmp(m.rxns, urea));
                F = fastingTargetFVA(m, {urea});  r.ureaFVA = [F.fvaMin F.fvaMax];
            end
            save(f, 'r');
        end
        fprintf('  urine creatinine %.1f mg/dL: feasibility %d, BMR %.0f kcal, urinary urea QP-FBA %.2f, FVA [%.2f %.2f]\n', ...
            r.uc, r.feasStat, r.kcal, r.ureaQP, r.ureaFVA);
    end
end
