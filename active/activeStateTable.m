function activeStateTable(res, runMETs, texFile)
% LaTeX table of activeWBM.m (Tables/activeState.tex): the model at rest, what the run
% needs, and for each scenario which limits were changed, the capacity for extra muscle
% work and the energy expenditure and RQ with the run imposed.
%
% INPUT: res, runMETs as saved by activeWBM.m in results/active/activeResults.mat
%
% Anna Sheehy 2026.

sexes = {'male', 'female'};
% scenario and the limits changed for the run
scen = {'diet', 'None'; 'stores', 'Stores released'; 'o2', 'Whole-body O$_2$ raised'; ...
        'o2stores', 'Whole-body O$_2$ raised, stores released'; 'o2muscle', 'Whole-body and muscle O$_2$ raised'; ...
        'o2musclestores', 'Whole-body and muscle O$_2$ raised, stores released'};
num  = @(x, f) ifelseStr(isfinite(x), sprintf(f, x), '--');

tf = fopen(texFile, 'w');
fprintf(tf, '\\begin{table}[H]\n\\centering\n\\small\n\\setlength{\\tabcolsep}{4pt}\n');
fprintf(tf, '\\begin{tabular}{@{}>{\\raggedright\\arraybackslash}p{4.6cm}rrrrrr@{}}\n\\toprule\n');
fprintf(tf, [' & \\multicolumn{2}{c}{\\textbf{Extra muscle work}} & ' ...
    '\\multicolumn{2}{c}{\\textbf{Energy expenditure}} & \\multicolumn{2}{c}{\\textbf{RQ}} \\\\\n']);
fprintf(tf, ' & \\multicolumn{2}{c}{(kcal/day)} & \\multicolumn{2}{c}{(kcal/day)} & & \\\\\n');
fprintf(tf, '\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}\n');
fprintf(tf, [' & \\textbf{Harvey} & \\textbf{Harvetta} & ' ...
    '\\textbf{Harvey} & \\textbf{Harvetta} & \\textbf{Harvey} & \\textbf{Harvetta} \\\\\n\\midrule\n']);

% model at rest and what the run needs
rq = @(f) f(2) / -f(1);                                   % fluxes(1:2) = O2 uptake, CO2 output
fprintf(tf, 'Model at rest & & & %.0f & %.0f & %.2f & %.2f \\\\\n', ...
    res.male.rest.kcal, res.female.rest.kcal, rq(res.male.rest.fluxes), rq(res.female.rest.fluxes));
fprintf(tf, 'Run, expected & %.0f & %.0f & %.0f & %.0f & & \\\\\n', ...
    res.male.boutKcal, res.female.boutKcal, res.male.rest.kcal + res.male.boutKcal, ...
    res.female.rest.kcal + res.female.boutKcal);
fprintf(tf, '\\midrule\n\\multicolumn{7}{@{}l}{\\textit{Model with the run; limits changed for the run:}} \\\\\n');

% scenarios
for r = 1:size(scen, 1)
    v = nan(2, 3);                                         % capacity, EE, RQ per sex
    for k = 1:2
        x = res.(sexes{k}).(scen{r, 1});
        v(k, 1) = x.capacityKcal;
        if ~isempty(x.fluxes), v(k, 2:3) = [x.kcal, rq(x.fluxes)]; end
    end
    fprintf(tf, '\\quad %s & %s & %s & %s & %s & %s & %s \\\\\n', scen{r, 2}, ...
        num(v(1, 1), '%.0f'), num(v(2, 1), '%.0f'), ...
        num(v(1, 2), '%.0f'), num(v(2, 2), '%.0f'), num(v(1, 3), '%.2f'), num(v(2, 3), '%.2f'));
end
fprintf(tf, '\\bottomrule\n\\end{tabular}\n');
fprintf(tf, ['\\caption[Capacity for extra muscle work and response to one hour of running]{\\textbf{Capacity for extra muscle work and response to one hour of running in the fed models.} ' ...
    'Run: one hour at %.1f METs, imposed as extra muscle ATP hydrolysis. ' ...
    'Extra muscle work: for the run, its cost; for the models, the largest extra muscle ATP hydrolysis they can support under BMR conditions (capacity, LP). ' ...
    'Energy expenditure and RQ: from the basal flux distribution (QP-FBA); --, the run exceeds the capacity. ' ...
    'Whole-body and muscle O$_2$ raised: the whole-body or muscle O$_2$ uptake limit was raised by the extra O$_2$ of the run (%.0f~mmol for Harvey and %.0f~mmol for Harvetta). ' ...
    'Stores released: hepatic and muscle glycogen and adipose triacylglycerol could be released.}\n'], ...
    runMETs, res.male.extraO2, res.female.extraO2);
fprintf(tf, '\\label{tab:active}\n\\end{table}\n');
fclose(tf);
end

function s = ifelseStr(c, a, b)
if c, s = a; else, s = b; end
end
