function oatEffectsFigure(long, rank, figPath, thesisFig, exclude)
% Figure of metabolomicsOAT.m (metabG_oatHeatmap): effect of each personalisation input on
% the feasible range of the 11 read-outs, Harvey (left) and Harvetta (right). Inputs with an
% effect of at least 0.01 on any read-out get their own row; the other metabolites share one
% row showing their largest effect. Inputs with an infeasible level and the height-only
% check are shown last. Cells give the effect (shift in feasible range / control range),
% labelled to one decimal where it is at least 0.05.
%
% INPUT: long     oatEffects_long.csv (sex, input, kind, readout, effect, ...)
%        rank     oatRanking.csv (rows ordered by score)
%        exclude  (optional) inputs reported separately, e.g. {'Urine creatinine'}, whose
%                 effect on urinary urea is too large for a shared colour scale
%
% Anna Sheehy 2026.

sexes = {'male', 'female'};  sexLab = {'Harvey', 'Harvetta'};
long.sex = string(long.sex);  long.input = string(long.input);  long.readout = string(long.readout);
readouts = unique(long.readout, 'stable');  nR = numel(readouts);
if nargin < 5, exclude = {}; end
rank = rank(~ismember(string(rank.input), string(exclude)), :);
inputs = string(rank.input);  kinds = string(rank.kind);
E = nan(numel(inputs), nR, 2);
for k = 1:2
    for i = 1:numel(inputs)
        for r = 1:nR
            v = long.effect(long.sex == sexes{k} & long.input == inputs(i) & long.readout == readouts(r));
            if ~isempty(v), E(i, r, k) = v(1); end
        end
    end
end
infeas = rank.infeasibleLevelMale | rank.infeasibleLevelFemale;
isHeight = inputs == "Height";
big = max(max(E, [], 3), [], 2) >= 0.01 & ~infeas & ~isHeight;
small = ~big & ~infeas & ~isHeight;
rows = {};  M = zeros(0, nR, 2);
for i = find(big)'
    rows{end+1, 1} = char(inputs(i)); M(end+1, :, :) = E(i, :, :); %#ok<AGROW>
end
if any(small)
    lab = sprintf('%d other metabolites (largest)', sum(small & kinds == "metabolite"));
    if any(small & kinds == "physiology"), lab = sprintf('%d other inputs (largest)', sum(small)); end
    rows{end+1, 1} = lab;  M(end+1, :, :) = max(E(small, :, :), [], 1);
end
for i = find(infeas)'                                    % height alone (noise check) is not shown
    rows{end+1, 1} = [char(inputs(i)) ' (low level infeasible)'];  M(end+1, :, :) = NaN; %#ok<AGROW>
end
rows = strrep(rows, 'Height (surface-area cardiac output)', 'Height, cardiac output scaled');
nRow = numel(rows);

cmax = max(M(:), [], 'omitnan');
cmap = interp1([0 1], [1 1 1; 0.03 0.19 0.42], linspace(0, 1, 256));
h = figure('Units', 'centimeters', 'Position', [2 2 16 0.5 * nRow + 5], 'Color', 'w');
tl = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for k = 1:2
    ax = nexttile;
    imagesc(M(:, :, k), 'AlphaData', ~isnan(M(:, :, k)));
    set(ax, 'Color', [0.85 0.85 0.85]);                   % infeasible rows grey
    colormap(ax, cmap);  caxis([0 cmax]);
    for i = 1:nRow
        for r = 1:nR
            v = M(i, r, k);
            if ~isnan(v) && v >= 0.05
                text(r, i, sprintf('%.1f', v), 'HorizontalAlignment', 'center', 'FontSize', 7, ...
                    'Color', [1 1 1] * (v > 0.55 * cmax));
            end
        end
    end
    set(ax, 'YTick', 1:nRow, 'YTickLabel', rows, 'XTick', 1:nR, 'XTickLabel', readouts, ...
        'XTickLabelRotation', 45, 'FontSize', 8, 'TickDir', 'out', 'TickLength', [0 0], 'Box', 'on');
    if k == 2, set(ax, 'YTickLabel', []); end
    title(sexLab{k}, 'FontWeight', 'normal');
end
cb = colorbar;  cb.Layout.Tile = 'east';
cb.Title.String = 'Effect';  cb.FontSize = 8;
exportgraphics(h, fullfile(figPath, 'metabG_oatHeatmap.png'), 'Resolution', 300);
exportgraphics(h, fullfile(figPath, 'metabG_oatHeatmap.pdf'), 'ContentType', 'vector');
if exist(thesisFig, 'dir'), copyfile(fullfile(figPath, 'metabG_oatHeatmap.png'), thesisFig); end
end
