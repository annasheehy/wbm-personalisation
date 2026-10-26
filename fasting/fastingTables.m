%% fastingTables.m
% Thesis Ch5, fasted state: two LaTeX tables from the existing caches of fastingWBM.m /
% fastingFigures.m (no new solves).
%   Tables/fastingWholeBody.tex   fed vs fasted BMR (L2), O2 uptake, CO2 output, RQ
%   Tables/fastingBenchmarks.tex  per benchmark and sex: literature range, basal (L2) flux,
%                                 feasible range (FVA) and whether each lies in the range
% Values in the benchmark table are oriented so that a positive number is
% the direction described in the literature (e.g. uptake shown as a positive uptake).
%
% Anna Sheehy 2026.

P = persWBMpaths();
resPath   = fullfile(P.results, 'fasting');
tabPath   = P.tables;
if ~exist(tabPath, 'dir'), mkdir(tabPath); end
sexes = {'male', 'female'};  sexLab = {'Harvey', 'Harvetta'};

load(fullfile(resPath, 'ckpt_figData.mat'), 'D');
load(fullfile(resPath, 'ckpt_fvaTargets.mat'), 'fvaM', 'fvaF');
FVA = struct('male', fvaM, 'female', fvaF);
Mdl = load(fullfile(resPath, 'ckpt_feasibleModels.mat'), 'maleFed', 'femaleFed', 'maleFast', 'femaleFast');
L2  = load(fullfile(resPath, 'ckpt_fluxL2.mat'), 'fluxFedM', 'fluxFedF', 'fluxFastM', 'fluxFastF');

%% Table 1: whole-body measures ---------------------------------------------------------
bmrOf = @(m, v) full(sum(([max(-m.S, 0), max(m.S, 0)] * [max(v, 0); max(-v, 0)]) .* contains(m.mets, '_atp['))) ...
    / 1000 * 64 * 0.239006;
W = struct();
pairs = {'male', 'maleFed', 'maleFast', 'fluxFedM', 'fluxFastM'; 'female', 'femaleFed', 'femaleFast', 'fluxFedF', 'fluxFastF'};
for k = 1:2
    s = pairs{k, 1};  d = D.(s);
    vFed = full(L2.(pairs{k, 4}));  vFast = full(L2.(pairs{k, 5}));
    W.(s).bmr = [bmrOf(Mdl.(pairs{k, 2}), vFed(1:numel(d.rxns))), bmrOf(Mdl.(pairs{k, 3}), vFast(1:numel(d.rxns)))];
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
Tt = fastingLiteratureTargets([], 'male');
rows = find(Tt.fluxComparable & ~isnan(Tt.litLow) & strlength(Tt.rxn) > 0);
B = struct();
for k = 1:2
    s = sexes{k};  T = fastingLiteratureTargets([], s);  d = D.(s);  f = FVA.(s);
    for n = 1:numel(rows)
        r = rows(n);
        i = find(strcmp(d.rxns, char(T.rxn(r))), 1);  j = find(f.rxn == T.rxn(r), 1);
        sgn = 1;  if strcmp(T.expectDir(r), 'uptake'), sgn = -1; end
        if ~d.posIsRelease(i), sgn = -sgn; end
        fv = sort(sgn * [f.fvaMin(j) f.fvaMax(j)]);
        B(n).label = char(T.label(r));  B(n).rxn = char(T.rxn(r));
        B(n).noSource = strlength(T.source(r)) == 0;  B(n).source = char(T.source(r));
        B(n).(s).lit = [T.litLow(r) T.litHigh(r)];
        B(n).(s).basal = sgn * d.vFast(i);
        B(n).(s).fva = fv;
        B(n).(s).verdict = plainVerdict(sgn * d.vFast(i), fv, [T.litLow(r) T.litHigh(r)]);
        [B(n).(s).symBasal, B(n).(s).symFeasible] = rangeSymbols(sgn * d.vFast(i), fv, [T.litLow(r) T.litHigh(r)]);
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
    'Measured: range reported in humans \\cite{%s}. All values are in the direction named, so an uptake is a positive uptake. ' ...
    'QP-FBA: flux in the basal flux distribution. FVA range: smallest and largest possible flux. ' ...
    'Within measured range: $\\checkmark$ within it, $\\downarrow$ below it, $\\uparrow$ above it; for FVA, $\\times$ means no possible flux lies in it. ' ...
    'Imposed: the measured range is a bound of the model, so the model matches it by construction.%s}\n'], ...
    strjoin(src, ','), repmat(' $^{\\ddagger}$No source assigned yet.', 1, any([B.noSource])));
fprintf(tf, '\\label{tab:fastBenchmarks}\n\\end{table}\n');
fclose(tf);

fprintf('Whole body: male BMR %.0f -> %.0f, female %.0f -> %.0f kcal/day\n', W.male.bmr, W.female.bmr);
for n = 1:numel(B), fprintf('%-36s M: %-22s F: %s\n', B(n).label, B(n).male.verdict, B(n).female.verdict); end

%% helpers -------------------------------------------------------------------------------
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
