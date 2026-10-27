function oatL2Figure(T, figPath, thesisFig)
% Figure of metabolomicsOAT.m (metabH_oatL2): for each personalisation input, the number of
% reactions whose basal flux changed between its low and high value (left) and the change in
% BMR (right), Harvey and Harvetta. Physiological inputs and metabolites that changed BMR by
% at least 5 kcal/day in either model get their own row, ordered by the mean absolute BMR
% change; the other metabolites share one row showing the range of their values. Inputs with
% a level without solution are shown last; the height-only check is not shown.
%
% INPUT: T  oatL2Ranking.csv
%
% Anna Sheehy 2026.

col = [0.13 0.40 0.67; 0.80 0.25 0.33];  mk = {'o', 's'};
bad = T.noSolutionMale | T.noSolutionFemale;
isHeight = strcmp(T.input, 'Height');
N = [T.nChangedMale, T.nChangedFemale];  D = [T.dBMRMale, T.dBMRFemale];
small = strcmp(T.kind, 'metabolite') & ~bad & max(abs(D), [], 2) < 5;
own = find(~small & ~bad & ~isHeight);
[~, o] = sort(mean(abs(D(own, :)), 2), 'descend');  own = own(o);
labels = strrep(T.input(own), 'Height (surface-area cardiac output)', 'Height, cardiac output scaled');
rowN = N(own, :);  rowD = D(own, :);  isRange = false(numel(own), 1);
if any(small)                                             % one row: range over the small metabolites
    labels{end+1} = sprintf('%d other metabolites (range)', sum(small));
    rowN(end+1, :) = NaN;  rowD(end+1, :) = NaN;  isRange(end+1) = true;
end
for i = find(bad)'                                       % height alone (noise check) is not shown
    labels{end+1} = [T.input{i} ' (low level infeasible)'];  rowN(end+1, :) = NaN;  rowD(end+1, :) = NaN;  isRange(end+1) = false; %#ok<AGROW>
end
n = numel(labels);  r0 = find(isRange);

h = figure('Units', 'centimeters', 'Position', [2 2 16 0.45 * n + 3], 'Color', 'w');
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for panel = 1:2
    nexttile; hold on
    if panel == 1, X = rowN; G = N(small, :); else, X = rowD; G = D(small, :); end
    for k = 1:2
        y = (1:n)' + (k - 1.5) * 0.2;
        pos = ~isnan(X(:, k)) & (panel == 2 | X(:, k) > 0);
        p(k) = plot(X(pos, k), y(pos), mk{k}, 'MarkerFaceColor', col(k, :), 'MarkerEdgeColor', 'none', 'MarkerSize', 4); %#ok<AGROW>
        if ~isempty(r0)
            plot([min(G(:, k)) max(G(:, k))], y(r0) * [1 1], '-', 'Color', col(k, :), 'LineWidth', 2);
        end
    end
    if panel == 1
        for j = find(all(X == 0, 2))', text(1.15, j, '0', 'FontSize', 8); end   % no reaction changed
        set(gca, 'XScale', 'log', 'XLim', [1 1e5], 'XTick', [1 100 1e4], ...
            'XTickLabel', {'1', '100', '10,000'}, 'YTickLabel', labels);
        xlabel('Reactions with a changed flux');
    else
        xline(0, 'k-');
        set(gca, 'YTickLabel', []);
        xlabel('BMR change, high minus low (kcal/day)');
    end
    hold off
    set(gca, 'YDir', 'reverse', 'YTick', 1:n, 'FontSize', 8, 'TickDir', 'out', 'Box', 'off', 'YLim', [0.4 n + 0.6]);
end
lg = legend(p, {'Harvey', 'Harvetta'}, 'Orientation', 'horizontal', 'Box', 'off');
lg.Layout.Tile = 'south';
exportgraphics(h, fullfile(figPath, 'metabH_oatL2.png'), 'Resolution', 300);
exportgraphics(h, fullfile(figPath, 'metabH_oatL2.pdf'), 'ContentType', 'vector');
if exist(thesisFig, 'dir'), copyfile(fullfile(figPath, 'metabH_oatL2.png'), thesisFig); end
end
