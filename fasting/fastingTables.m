%% fastingTables.m
% Thesis Ch5, fasted state: two LaTeX tables from the models and L2 fluxes of fastingWBM.m.
%   Tables/fastingWholeBody.tex   fed vs fasted BMR (L2), O2 uptake, CO2 output, RQ
%   Tables/fastingBenchmarks.tex  per benchmark and sex: literature range, basal (L2) flux,
%                                 feasible range (FVA) and whether each lies in the range
% Values in the benchmark table are oriented so that a positive number is the direction
% described in the literature (e.g. uptake shown as a positive uptake). Bounds, fluxes and
% the FVA of the benchmark reactions are cached (ckpt_figData.mat, ckpt_fvaTargets.mat), so
% only the first run solves anything (FVA, ~1 min per target).
%
% Anna Sheehy 2026.

P = persWBMpaths();
resPath   = fullfile(P.results, 'fasting');
tabPath   = P.tables;
if ~exist(tabPath, 'dir'), mkdir(tabPath); end
sexes = {'male', 'female'};  sexLab = {'Harvey', 'Harvetta'};

figData = fullfile(resPath, 'ckpt_figData.mat');
fvaData = fullfile(resPath, 'ckpt_fvaTargets.mat');
if exist(figData, 'file'), load(figData, 'D'); else, D = benchmarkFluxData(resPath); save(figData, 'D'); end
% FVA of every benchmark reaction (or signed sum) under BMR conditions; only missing ones are solved
Tt = fastingLiteratureTargets([], 'male');
rl = unique(cellstr(Tt.rxn(strlength(Tt.rxn) > 0)), 'stable');
fvaM = [];  fvaF = [];
if exist(fvaData, 'file'), load(fvaData, 'fvaM', 'fvaF'); end
todo = rl;  if ~isempty(fvaM), todo = rl(~ismember(rl, cellstr(fvaM.rxn))); end
if ~isempty(todo)
    changeCobraSolver('gurobi', 'LP', 0);
    mdl = load(fullfile(resPath, 'ckpt_feasibleModels.mat'), 'maleFast', 'femaleFast');
    fvaM = [fvaM; fastingTargetFVA(mdl.maleFast, todo)];
    fvaF = [fvaF; fastingTargetFVA(mdl.femaleFast, todo)];
    save(fvaData, 'fvaM', 'fvaF');
end
FVA = struct('male', fvaM, 'female', fvaF);
Mdl = load(fullfile(resPath, 'ckpt_feasibleModels.mat'), 'maleFed', 'femaleFed', 'maleFast', 'femaleFast');
L2  = load(fullfile(resPath, 'ckpt_fluxL2.mat'), 'fluxFedM', 'fluxFedF', 'fluxFastM', 'fluxFastF');

%% Table 1: whole-body measures ---------------------------------------------------------
W = struct();
pairs = {'male', 'maleFed', 'maleFast', 'fluxFedM', 'fluxFastM'; 'female', 'femaleFed', 'femaleFast', 'fluxFedF', 'fluxFastF'};
for k = 1:2
    s = pairs{k, 1};  d = D.(s);
    vFed = full(L2.(pairs{k, 4}));  vFast = full(L2.(pairs{k, 5}));
    W.(s).bmr = [bmrFromFlux(Mdl.(pairs{k, 2}), vFed), bmrFromFlux(Mdl.(pairs{k, 3}), vFast)];
    W.(s).o2  = -[d.vFed(d.iO2), d.vFast(d.iO2)];
    W.(s).co2 = [d.vFed(d.iCO2), d.vFast(d.iCO2)];
    W.(s).rq  = W.(s).co2 ./ W.(s).o2;
    W.(s).o2max = -d.lbFast(d.iO2);  W.(s).co2min = d.lbFast(d.iCO2);
end
tf = fopen(fullfile(tabPath, 'fastingWholeBody.tex'), 'w');
fprintf(tf, '\\begin{table}[H]\n\\centering\n\\small\n');
fprintf(tf, '\\begin{tabular}{@{}lrrrr@{}}\n\\toprule\n');
fprintf(tf, ' & \\multicolumn{2}{c}{\\textbf{%s}} & \\multicolumn{2}{c}{\\textbf{%s}} \\\\\n\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\n', sexLab{:});
fprintf(tf, ' & \\textbf{Fed} & \\textbf{Fasted} & \\textbf{Fed} & \\textbf{Fasted} \\\\\n\\midrule\n');
rowsW = {'BMR (kcal per day)', 'bmr', '%.0f'; 'O$_2$ uptake (mmol per day)', 'o2', '%.0f'; ...
         'CO$_2$ output (mmol per day)', 'co2', '%.0f'; 'Respiratory quotient', 'rq', '%.2f'};
for r = 1:size(rowsW, 1)
    f = rowsW{r, 2};  fm = rowsW{r, 3};
    fprintf(tf, ['%s & ' fm ' & ' fm ' & ' fm ' & ' fm ' \\\\\n'], rowsW{r, 1}, W.male.(f)(1), W.male.(f)(2), ...
        W.female.(f)(1), W.female.(f)(2));
end
fprintf(tf, '\\bottomrule\n\\end{tabular}\n');
fprintf(tf, ['\\caption[BMR and gas exchange of the fed and fasted models]{\\textbf{BMR and gas exchange of the fed and fasted models.} ' ...
    'All values from the basal flux distribution (QP-FBA). The O$_2$ uptake limit is %.0f~mmol per day and the lower limit on CO$_2$ output %.0f~mmol per day in both models.}\n'], ...
    W.male.o2max, W.male.co2min);
fprintf(tf, '\\label{tab:fastWholeBody}\n\\end{table}\n');
fclose(tf);

%% Table 2: benchmarks ------------------------------------------------------------------
rows = find(Tt.fluxComparable & ~isnan(Tt.litLow) & strlength(Tt.rxn) > 0);
B = struct();
for k = 1:2
    s = sexes{k};  T = fastingLiteratureTargets([], s);  d = D.(s);  f = FVA.(s);
    for n = 1:numel(rows)
        r = rows(n);
        [i, c] = netExchange(d.rxns, T.rxn(r));  j = find(f.rxn == T.rxn(r), 1);
        sgn = 1;  if strcmp(T.expectDir(r), 'uptake'), sgn = -1; end
        if isscalar(i) && ~d.posIsRelease(i), sgn = -sgn; end     % single reaction written the other way round
        v = sgn * (c' * d.vFast(i));                              % a signed sum is already net release
        fv = sort(sgn * [f.fvaMin(j) f.fvaMax(j)]);
        B(n).label = char(T.label(r));  B(n).rxn = char(T.rxn(r));
        B(n).noSource = strlength(T.source(r)) == 0;  B(n).source = char(T.source(r));
        B(n).(s).lit = [T.litLow(r) T.litHigh(r)];
        B(n).(s).basal = v;
        B(n).(s).fva = fv;
        B(n).(s).verdict = plainVerdict(v, fv, [T.litLow(r) T.litHigh(r)]);
        [B(n).(s).symBasal, B(n).(s).symFeasible] = rangeSymbols(v, fv, [T.litLow(r) T.litHigh(r)]);
    end
end
tf = fopen(fullfile(tabPath, 'fastingBenchmarks.tex'), 'w');
fprintf(tf, '\\begin{table}[H]\n\\centering\n\\small\n\\setlength{\\tabcolsep}{4pt}\n');
fprintf(tf, '\\begin{tabular}{@{}lrrrcc@{}}\n\\toprule\n');
fprintf(tf, ' & & \\multicolumn{2}{c}{\\textbf{Model}} & \\multicolumn{2}{c}{\\textbf{Within measured range}} \\\\\n\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}\n');
fprintf(tf, '\\textbf{Flux (mmol per day)} & \\textbf{Measured} & \\textbf{QP-FBA} & \\textbf{FVA range} & \\textbf{QP-FBA} & \\textbf{FVA} \\\\\n\\midrule\n');
for k = 1:2
    s = sexes{k};
    if k > 1, fprintf(tf, '\\midrule\n'); end
    fprintf(tf, '\\multicolumn{6}{@{}l}{\\textit{%s}} \\\\\n', sexLab{k});
    for n = 1:numel(B)
        v = B(n).(s);
        lab = B(n).label;
        if strcmp(v.verdict, 'set by a model bound'), lab = [lab ' (imposed)']; end
        if B(n).noSource, lab = [lab '$^{\ddagger}$']; end
        fprintf(tf, '%s & %s--%s & %s & %s to %s & %s & %s \\\\\n', lab, num(v.lit(1)), num(v.lit(2)), ...
            num(v.basal), num(v.fva(1)), num(v.fva(2)), v.symBasal, v.symFeasible);
    end
end
fprintf(tf, '\\bottomrule\n\\end{tabular}\n');
src = unique({B(~[B.noSource]).source}, 'stable');
fprintf(tf, ['\\caption[Fasted models compared with fluxes measured in humans]{\\textbf{Organ exchange fluxes of the fasted models compared with the same fluxes measured in humans after an overnight to 24~h fast.} ' ...
    'Measured: range reported in humans \\cite{%s}. For red blood cell glucose uptake and muscle and hepatic alanine exchange, the range is a single reported rate $\\pm 20\\%%$; for brain glucose uptake, the lower end is set 20\\%% below the reported rate. All values are in the direction named, so an uptake is a positive uptake and a negative value means the opposite direction. ' ...
    'Endogenous glucose production (liver and kidney), renal glucose release and hepatic alanine uptake are net exchanges over all blood compartments. ' ...
    'QP-FBA: flux in the basal flux distribution. FVA range: smallest and largest possible flux. ' ...
    'Within measured range: $\\checkmark$ within it, $\\downarrow$ below it, $\\uparrow$ above it; for FVA, $\\times$ means no possible flux lies in it. ' ...
    'Imposed: the measured range is a bound of the model, so the model matches it by construction.%s}\n'], ...
    strjoin(src, ','), repmat(' $^{\\ddagger}$No source assigned yet.', 1, any([B.noSource])));
fprintf(tf, '\\label{tab:fastBenchmarks}\n\\end{table}\n');
fclose(tf);

fprintf('Whole body: male BMR %.0f -> %.0f, female %.0f -> %.0f kcal/day\n', W.male.bmr, W.female.bmr);
for n = 1:numel(B), fprintf('%-36s M: %-22s F: %s\n', B(n).label, B(n).male.verdict, B(n).female.verdict); end

%% helpers -------------------------------------------------------------------------------
function D = benchmarkFluxData(resPath)
% Reaction ids, bounds and L2 fluxes of the fed and fasted models, read once from the
% checkpoints of fastingWBM.m. posIsRelease records, per reaction and from the
% stoichiometry, whether a positive flux is release out of the organ (into [bc], [csf],
% [u] or [a]).
load(fullfile(resPath, 'ckpt_feasibleModels.mat'), 'maleFed', 'femaleFed', 'maleFast', 'femaleFast');
load(fullfile(resPath, 'ckpt_fluxL2.mat'), 'fluxFedM', 'fluxFedF', 'fluxFastM', 'fluxFastF');
sets = {'male', maleFed, maleFast, fluxFedM, fluxFastM; 'female', femaleFed, femaleFast, fluxFedF, fluxFastF};
for s = 1:2
    mFed = sets{s, 2};  mFast = sets{s, 3};
    d = struct('rxns', {mFast.rxns}, 'lbFed', mFed.lb, 'ubFed', mFed.ub, 'lbFast', mFast.lb, 'ubFast', mFast.ub, ...
        'vFed', full(sets{s, 4}), 'vFast', full(sets{s, 5}));
    d.iO2  = find(strcmp(d.rxns, 'EX_o2[a]'), 1);
    d.iCO2 = find(strcmp(d.rxns, 'EX_co2[a]'), 1);
    d.posIsRelease = true(numel(d.rxns), 1);
    for i = find(~cellfun(@isempty, regexp(d.rxns, '_EX_|^EX_', 'once')))'
        mr  = find(mFast.S(:, i) ~= 0);
        out = mr(~cellfun(@isempty, regexp(mFast.mets(mr), '\[(bc|csf|u|a)\]$', 'once')));
        if isempty(out), continue; end
        d.posIsRelease(i) = full(mFast.S(out(1), i)) > 0 || (numel(mr) == 1 && full(mFast.S(mr, i)) < 0);
    end
    D.(sets{s, 1}) = d;
end
end

function v = plainVerdict(basal, fva, lit)
tol = 1e-6 * max(1, max(abs(lit)));
if fva(1) >= lit(1) - tol && fva(2) <= lit(2) + tol
    v = 'set by a model bound';
elseif basal >= lit(1) - tol && basal <= lit(2) + tol
    v = 'within range';
elseif fva(2) >= lit(1) && fva(1) <= lit(2)
    v = 'possible, not chosen';
else
    v = 'not possible';
end
end

function [b, f] = rangeSymbols(basal, fva, lit)
% basal flux vs the literature range (tick / below / above); whether any feasible flux lies in it
tol = 1e-6 * max(1, max(abs(lit)));
if basal < lit(1) - tol,     b = '$\downarrow$';
elseif basal > lit(2) + tol, b = '$\uparrow$';
else,                        b = '$\checkmark$';
end
if fva(2) >= lit(1) - tol && fva(1) <= lit(2) + tol, f = '$\checkmark$'; else, f = '$\times$'; end
end

function s = num(x)
if ~isfinite(x), s = '--'; elseif abs(x) >= 10, s = sprintf('%.0f', x); elseif abs(x) >= 1, s = sprintf('%.1f', x); else, s = sprintf('%.2g', x); end
if x < 0, s = strrep(s, '-', '$-$'); end           % minus sign, not a hyphen
end
