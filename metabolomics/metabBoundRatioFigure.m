function metabBoundRatioFigure(R, metNames, w, figPath, thesisFig)
% Figure of metabolomicsPersonalisation.m (metabB_boundRatio): organ uptake limits after
% personalisation relative to the default limits, in Harvey (left) and Harvetta (right).
% One box per metabolite over the synthetic profiles; metabolites in the same order in both
% panels (median ratio over both sexes, largest at the top).
%
% INPUT: R.(sex)   nProfiles x nMet, median over organs of personalised / default uptake limit
%        metNames  display names, same order as the columns of R
%        w         measured concentration allowance used (for the console only)
%
% Anna Sheehy 2026.

sexes = {'male', 'female'};  sexLab = {'Harvey', 'Harvetta'};
col = struct('male', [0.13 0.40 0.67], 'female', [0.80 0.25 0.33]);
nMet = numel(metNames);
[~, ord] = sort(median([R.male; R.female], 1, 'omitnan'));        % smallest at the bottom
h = figure('Units', 'centimeters', 'Position', [2 2 16 0.33 * nMet + 2.5], 'Color', 'w');
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
xl = [min([R.male(:); R.female(:)]), max([R.male(:); R.female(:)])] .* [0.8 1.25];
for k = 1:2
    s = sexes{k};  X = R.(s)(:, ord);
    g = repmat(1:nMet, size(X, 1), 1);
    nexttile;
    boxchart(g(:), X(:), 'Orientation', 'horizontal', 'BoxFaceColor', col.(s), 'MarkerColor', col.(s), ...
        'MarkerStyle', '.', 'BoxWidth', 0.6);
    set(gca, 'XScale', 'log', 'YTick', 1:nMet, 'YTickLabel', metNames(ord), 'FontSize', 8, ...
        'TickDir', 'out', 'Box', 'off', 'YLim', [0.4 nMet + 0.6], 'XLim', xl);
    tk = [0.1 0.2 0.5 1 2 5 10];  tk = tk(tk >= xl(1) & tk <= xl(2));
    set(gca, 'XTick', tk, 'XTickLabel', compose('%g', tk));
    if k == 2, set(gca, 'YTickLabel', []); end
    xline(1, 'k--');
    xlabel('Personalised / default uptake limit');
    title(sexLab{k}, 'FontWeight', 'normal');
end
fprintf('metabBoundRatioFigure: allowance +/-%d%%, %d profiles per sex\n', round(100 * w), size(R.male, 1));
exportgraphics(h, fullfile(figPath, 'metabB_boundRatio.png'), 'Resolution', 300);
exportgraphics(h, fullfile(figPath, 'metabB_boundRatio.pdf'), 'ContentType', 'vector');
if exist(thesisFig, 'dir'), copyfile(fullfile(figPath, 'metabB_boundRatio.png'), thesisFig); end
end
