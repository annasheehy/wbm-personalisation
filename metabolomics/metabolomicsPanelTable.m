function metabolomicsPanelTable(panelCsv, texFile, nPerSex, defaultConc)
% LaTeX table of the blood metabolite panel of metabolomicsPersonalisation.m, written to
% texFile as one table in two parts with the same columns: metabolites whose source reports
% men and women together (one healthy range, used for both sexes), then metabolites
% reported by sex. Each row gives the healthy range (2.5th-97.5th percentile of the fitted
% log-normal) and the concentration from which the default uptake bounds of the control
% models are calculated. Sources are cited in the caption, once per source.
%
% INPUT: panelCsv    inputPanel.csv written by metabolomicsPersonalisation.m (Part 8)
%        nPerSex     synthetic profiles per sex
%        defaultConc toolbox default blood concentration (uM), used where HMDB has no value
%
% Anna Sheehy 2026.

P = readtable(panelCsv, 'TextType', 'char', 'Delimiter', ',');
P.sexSpecific = logical(P.sexSpecific);
P.hmdbBelowHealthy = logical(P.hmdbBelowHealthy);
P.toolboxDefault = logical(P.toolboxDefault);
G1 = P(~P.sexSpecific, :);  G2 = P(P.sexSpecific, :);
hdr = @(sub) sprintf([' & & \\multicolumn{2}{c}{\\textbf{Healthy range ($\\mu$M)}} & \\textbf{Upper limit in} \\\\\n' ...
    '\\cmidrule(lr){3-4}\n\\textbf{Metabolite} & \\textbf{VMH ID} & %s & \\textbf{default models ($\\mu$M)} \\\\\n\\midrule\n'], sub);

tf = fopen(texFile, 'w');
fprintf(tf, '\\begin{table}[H]\n\\centering\n\\small\n');
fprintf(tf, '\\begin{tabular}{@{}llrrr@{}}\n\\toprule\n');
fprintf(tf, '%s', hdr('\multicolumn{2}{c}{\textbf{Men and women}}'));
for q = 1:height(G1)
    fprintf(tf, '%s & \\texttt{%s} & \\multicolumn{2}{c}{%s--%s} & %s \\\\\n', G1.name{q}, ...
        strrep(G1.vmh{q}, '_', '\_'), fmtNum(G1.p025Male(q)), fmtNum(G1.p975Male(q)), limitCell(G1(q, :), defaultConc));
end
if ~isempty(G2)
    fprintf(tf, '\\bottomrule\n\\addlinespace[10pt]\n\\toprule\n');
    fprintf(tf, '%s', hdr('\textbf{Men} & \textbf{Women}'));
    for q = 1:height(G2)
        fprintf(tf, '%s & \\texttt{%s} & %s--%s & %s--%s & %s \\\\\n', G2.name{q}, ...
            strrep(G2.vmh{q}, '_', '\_'), fmtNum(G2.p025Male(q)), fmtNum(G2.p975Male(q)), ...
            fmtNum(G2.p025Female(q)), fmtNum(G2.p975Female(q)), limitCell(G2(q, :), defaultConc));
    end
end
fprintf(tf, '\\bottomrule\n\\end{tabular}\n');
head = 'Healthy blood concentrations of the personalised metabolites';
fprintf(tf, '\\caption[%s]{\\textbf{%s.} Upper part: metabolites reported for men and women together%s; one range is used for both sexes. ', head, head, citeList(G1));
if ~isempty(G2)
    fprintf(tf, 'Lower part: metabolites reported for men and women separately%s. ', citeList(G2));
end
fprintf(tf, ['Healthy range: 2.5th to 97.5th percentile of the log-normal distribution fitted to the reported values; the %d synthetic profiles per sex were drawn from these distributions. ' ...
    'Upper limit in default models: the concentration from which \\texttt{physiologicalConstraintsHMDBbased} sets the default organ uptake limits, the highest healthy blood concentration in HMDB. ' ...
    '(default): no HMDB value, so the default of %g~$\\mu$M is used. Bold: below the upper end of the healthy range. VMH ID: Virtual Metabolic Human identifier.}\n'], nPerSex, defaultConc);
fprintf(tf, '\\label{tab:metabPanel}\n\\end{table}\n');
fclose(tf);
end

function s = citeList(G)
% ' \cite{A}' for one source; ', X and Y \cite{A}; Z \cite{B}' for several
src = unique(G.source, 'stable');
if numel(src) == 1, s = sprintf(' \\cite{%s}', src{1}); return, end
parts = cell(1, numel(src));
for k = 1:numel(src)
    nm = cellfun(@lowerName, G.name(strcmp(G.source, src{k})), 'UniformOutput', false);
    if numel(nm) > 1, nm = [strjoin(nm(1:end-1), ', ') ' and ' nm{end}]; else, nm = nm{1}; end
    parts{k} = sprintf('%s \\cite{%s}', nm, src{k});
end
s = [', ' strjoin(parts, '; ')];
end

function n = lowerName(n)
% metabolite name for running text: 'D-Glucose' -> 'D-glucose', 'Urea' -> 'urea'
if numel(n) > 2 && any(strcmp(n(1:2), {'D-', 'L-'})), n(3) = lower(n(3)); else, n(1) = lower(n(1)); end
end

function s = limitCell(R, defaultConc)
% concentration behind the default bounds; bold where it lies below the healthy range
if R.toolboxDefault
    v = sprintf('%g', defaultConc);
else
    v = fmtNum(R.hmdbUpperLimit);
end
if R.hmdbBelowHealthy, v = ['\textbf{' v '}']; end
if R.toolboxDefault, v = [v ' (default)']; end
s = v;
end

function s = fmtNum(x)
% Number for a LaTeX table: 0 decimals from 100, 1 from 10, else 2.
if ~isfinite(x), s = '--';
elseif abs(x) >= 100, s = sprintf('%.0f', x);
elseif abs(x) >= 10,  s = sprintf('%.1f', x);
else,                 s = sprintf('%.2f', x);
end
end
