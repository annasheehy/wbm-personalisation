%% fastingAdiposeFix.m
% Thesis Ch5, fasted state: tests constraints that make adipose tissue release the glycerol
% of lipolysis. In the fasted models of fastingWBM.m, adipose tissue exports stored
% triacylglycerol intact and takes up glycerol from the blood. Each variant adds direction
% constraints to the fasted models and recomputes the basal flux distribution (QP-FBA),
% BMR and the flux range (FVA) of the benchmark fluxes.
%   asBuilt   fasted models of fastingWBM.m
%   glyc      adipose glycerol exchange release only (no uptake)
%   glycTag   glyc + no adipose triacylglycerol release to blood
%   glycTagLiver  glycTag + hepatic glycerol exchange uptake only (no release)
%   frayn     adipose glycerol release of at least 0.8 times the postabsorptive rate of
%             Frayn, the lower bound written but commented out in
%             physiologicalConstraintsHMDBbased (its upper bound is 1.2 times that rate)
%   tagStore  adipose fat released only as triacylglycerol: the free fatty acid and acyl-CoA
%             stores opened by setFeedingFastingConstraints are closed, the triacylglycerol
%             store may release up to 1.2 times the Frayn rate (the adipose glycerol upper
%             bound), and glycTag applies
%   tagStoreLiver  tagStore + hepatic glycerol exchange uptake only (no release)
% OUTPUT: results/fasting/adiposeFix/<sex>_<variant>.mat (feasibility, BMR, fluxes, FVA)
%
% Anna Sheehy 2026.

P = persWBMpaths();
resPath = fullfile(P.results, 'fasting');
outPath = fullfile(resPath, 'adiposeFix');
if ~exist(outPath, 'dir'), mkdir(outPath); end
changeCobraSolver('gurobi', 'LP', 0);  changeCobraSolver('gurobi', 'QP', 0);

Mdl = load(fullfile(resPath, 'ckpt_feasibleModels.mat'), 'maleFast', 'femaleFast');
load(fullfile(resPath, 'ckpt_figData.mat'), 'D');
sexes = {'male', 'female'};
% variant, adipose glycerol lb, adipose TAG ub, liver glycerol ub (NaN = unchanged),
% fat stored only as TAG; 'frayn' = 0.8/1.2 of the adipose glycerol upper bound
variants = {'glycTag', 0, 0, NaN, false; 'glyc', 0, NaN, NaN, false; 'glycTagLiver', 0, 0, 0, false; ...
            'frayn', 'frayn', NaN, NaN, false; 'tagStore', 0, 0, NaN, true; 'tagStoreLiver', 0, 0, 0, true};
% adipose stores of free fatty acids and acyl-CoAs in setFeedingFastingConstraints
faStores = strcat('Adipocytes_sink_', {'hdca', 'tmndnc', 'lnlc', 'c226coa', 'doco13ecoa', 'lnlccoa', ...
    'lnlncacoa', 'lnlncgcoa', 'odecoa', 'pmtcoa', 'stcoa', 'tmndnccoa'}, '(c)');
show = {'Adipocytes_EX_glyc(e)_[bc]', 'Adipocytes_EX_tag_hs(e)_[bc]', 'Adipocytes_sink_tag_hs(c)', ...
        'Adipocytes_GLYK', 'Adipocytes_EX_hdca(e)_[bc]', 'Adipocytes_EX_ocdcea(e)_[bc]', ...
        'Liver_EX_glyc(e)_[bc]', 'Liver_EX_glc_D(e)_[bc]', 'Liver_EX_bhb(e)_[bc]', ...
        'Brain_EX_glc_D(e)_[csf]', 'Muscle_EX_ala_L(e)_[bc]'};

for k = 1:2
    s = sexes{k};  base = Mdl.([s 'Fast']);
    T = fastingLiteratureTargets([], s);
    bench = cellstr(T.rxn(T.fluxComparable & ~isnan(T.litLow) & strlength(T.rxn) > 0));
    gas = D.(s).rxns([D.(s).iO2, D.(s).iCO2]);
    for q = 1:size(variants, 1)
        f = fullfile(outPath, sprintf('%s_%s.mat', s, variants{q, 1}));
        if exist(f, 'file'), fprintf('%s %s: exists, skipped\n', s, variants{q, 1}); continue, end
        m = base;
        lbG = variants{q, 2};
        if strcmp(lbG, 'frayn'), lbG = m.ub(strcmp(m.rxns, 'Adipocytes_EX_glyc(e)_[bc]')) * 0.8 / 1.2; end
        if ~isnan(lbG), m = changeRxnBounds(m, 'Adipocytes_EX_glyc(e)_[bc]', lbG, 'l'); end
        if ~isnan(variants{q, 3}), m = changeRxnBounds(m, 'Adipocytes_EX_tag_hs(e)_[bc]', variants{q, 3}, 'u'); end
        if variants{q, 5}
            m = changeRxnBounds(m, faStores(ismember(faStores, m.rxns)), 0, 'l');
            m = changeRxnBounds(m, 'Adipocytes_sink_tag_hs(c)', -m.ub(strcmp(m.rxns, 'Adipocytes_EX_glyc(e)_[bc]')), 'l');
        end
        if ~isnan(variants{q, 4}), m = changeRxnBounds(m, 'Liver_EX_glyc(e)_[bc]', variants{q, 4}, 'u'); end
        sol = optimizeWBModel(m);
        r = struct('variant', variants{q, 1}, 'feasStat', sol.stat);
        fprintf('%s %s: feasibility stat %d\n', s, variants{q, 1}, sol.stat);
        if sol.stat == 1
            b = calculateModelBMR(m, 'norms', {'L2'});
            v = b.L2.v(1:numel(m.rxns));
            r.kcal = b.L2.kcal;  r.qpStat = b.L2.stat;
            r.show = show;  r.vShow = cellfun(@(x) v(strcmp(m.rxns, x)), show);
            r.gas = gas;  r.vGas = cellfun(@(x) v(strcmp(m.rxns, x)), gas);
            r.fva = fastingTargetFVA(m, unique([bench; show(1)], 'stable'));
            fprintf('%s %s: BMR %.0f kcal (stat %d), O2 %.0f, CO2 %.0f\n', s, variants{q, 1}, r.kcal, r.qpStat, r.vGas);
            for i = 1:numel(show), fprintf('   %-32s %10.1f\n', show{i}, r.vShow(i)); end
        end
        save(f, 'r');
    end
end
