function physiologyRangesTable(sourceFile, ipar, texFile)
% LaTeX table of the ranges used for the physiological inputs of metabolomicsOAT.m,
% with the default values of Harvey and Harvetta and the sources.
%
% Anna Sheehy 2026.

lims = healthyPhysiologyRanges(sourceFile, ipar);
S = readtable(sourceFile, 'TextType', 'char', 'Delimiter', ',');
rows = {'bodyWeight', 'Body weight (kg)', '%.1f'
        'Height', 'Height (cm)', '%.0f'
        'HeartRate', 'Resting heart rate (beats/min)', '%.0f'
        'Hematocrit', 'Haematocrit', '%.3f'
        'GlomerularFiltrationRate', 'Glomerular filtration rate (mL/min)', '%.0f'};
tf = fopen(texFile, 'w');
fprintf(tf, '\\begin{table}[H]\n\\centering\n\\small\n');
fprintf(tf, '\\begin{tabular}{@{}>{\\raggedright\\arraybackslash}p{3.4cm}rrrrl@{}}\n\\toprule\n');
fprintf(tf, ' & \\multicolumn{2}{c}{\\textbf{Men}} & \\multicolumn{2}{c}{\\textbf{Women}} & \\\\\n\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\n');
fprintf(tf, '\\textbf{Input} & \\textbf{Model default} & \\textbf{Range tested} & \\textbf{Model default} & \\textbf{Range tested} & \\textbf{Source} \\\\\n\\midrule\n');
for r = 1:size(rows, 1)
    f = rows{r, 1};  fm = rows{r, 3};
    src = S.source(strcmp(S.field, f));  src = src{1};
    if strcmp(src, 'TBD'), cit = '\todo{source}'; else, cit = sprintf('\\cite{%s}', src); end
    cells = cell(1, 2);
    for k = 1:2
        s = {'male', 'female'};  s = s{k};
        cells{k} = sprintf([fm ' & ' fm '--' fm], ipar.(s).(f), lims.(s).(f)(1), lims.(s).(f)(2));
    end
    fprintf(tf, '%s & %s & %s & %s \\\\\n', rows{r, 2}, cells{1}, cells{2}, cit);
end
fprintf(tf, '\\bottomrule\n\\end{tabular}\n');
fprintf(tf, '\\caption[Ranges tested for the physiological inputs]{\\textbf{Ranges tested for the physiological inputs.} Model default: value in Harvey (men) and Harvetta (women). Range tested: low and high values used for each input. Body weight: body mass index 18.5 to 24.9~kg\\,m$^{-2}$ at the reference height. Height: adults aged 20 to 39 years, mean $\\pm$ 1.96 SD, with the SD derived from the reported 5th and 95th percentiles. Glomerular filtration rate: estimated GFR of healthy adults aged 18 to 39 years, mean $\\pm$ 1.96 SD, converted from mL\\,min$^{-1}$ per 1.73~m$^2$ with the body surface area of the reference individual.}\n');
fprintf(tf, '\\label{tab:physRanges}\n\\end{table}\n');
fclose(tf);
end
