function groups = metabBoundMapFigure(map, metIDs, metNames, figPath, thesisFig)
% Figure of metabolomicsPersonalisation.m (metabA_boundMap): reactions constrained by each
% personalised blood metabolite, by type, in Harvey (left) and Harvetta (right).
% Metabolites with identical counts of every type in both models are shown as one bar,
% labelled with the first non-amino acid of the group (if any) and the number of amino acids.
%
% INPUT:  map       bound map of metabolomicsPersonalisation.m (ckpt_map.mat), map.(sex)
%         metIDs    VMH IDs of the panel; metNames = display names in the same order
% OUTPUT: groups    one row per bar: label, members (cell of names), counts per sex and type
%
% Anna Sheehy 2026.

sexes = {'male', 'female'};  sexLab = {'Harvey', 'Harvetta'};
typeOrder = {'organ uptake', 'kidney filtration', 'blood-brain barrier'};
typeLab = {'Organ uptake from blood', 'Kidney filtration', 'Blood to cerebrospinal fluid'};
nMet = numel(metIDs);
cnt = zeros(nMet, numel(typeOrder), 2);
for k = 1:2
    mp = map.(sexes{k});
    for q = 1:nMet
        for t = 1:numel(typeOrder)
            cnt(q, t, k) = numel(unique(mp.rxnIdx(strcmp(mp.met, metIDs{q}) & strcmp(mp.type, typeOrder{t}))));
        end
    end
end

% group metabolites with identical counts in both models
key = reshape(cnt, nMet, []);
[~, ~, g] = unique(key, 'rows', 'stable');
nG = max(g);
label = cell(nG, 1);  members = cell(nG, 1);  C = zeros(nG, numel(typeOrder), 2);
isAA = startsWith(metNames(:), 'L-') & ~strcmp(metNames(:), 'L-Carnitine') | strcmp(metNames(:), 'Glycine');
for j = 1:nG
    in = find(g == j);
    members{j} = metNames(in);
    C(j, :, :) = cnt(in(1), :, :);
    if numel(in) == 1
        label{j} = metNames{in};
    else
        nAA = sum(isAA(in));  other = metNames(in(~isAA(in)));
        lab = sprintf('%d amino acids', nAA);
        if ~isempty(other), lab = [strjoin(other(:)', ', ') ' and ' lab]; end
        label{j} = lab;
    end
end
[~, ord] = sort(sum(C(:, :, 1), 2) + sum(C(:, :, 2), 2));
label = label(ord);  members = members(ord);  C = C(ord, :, :);
groups = table(label, members, C(:, :, 1), C(:, :, 2), 'VariableNames', {'label', 'members', 'male', 'female'});

h = figure('Units', 'centimeters', 'Position', [2 2 16 0.55 * nG + 3], 'Color', 'w');
tl = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
xMax = ceil(max(sum(C, 2), [], 'all') / 5) * 5;
for k = 1:2
    nexttile;
    b = barh(C(:, :, k), 'stacked', 'BarWidth', 0.7);
    set(gca, 'YTick', 1:nG, 'YTickLabel', label, 'FontSize', 8, 'TickDir', 'out', 'Box', 'off', ...
        'XLim', [0 xMax], 'YLim', [0.4 nG + 0.6]);
    if k == 2, set(gca, 'YTickLabel', []); end
    xlabel('Reactions constrained');
    title(sexLab{k}, 'FontWeight', 'normal');
    if k == 2
        lg = legend(b, typeLab, 'Orientation', 'horizontal', 'Box', 'off');
        lg.Layout.Tile = 'south';
    end
end
exportgraphics(h, fullfile(figPath, 'metabA_boundMap.png'), 'Resolution', 300);
exportgraphics(h, fullfile(figPath, 'metabA_boundMap.pdf'), 'ContentType', 'vector');
if exist(thesisFig, 'dir'), copyfile(fullfile(figPath, 'metabA_boundMap.png'), thesisFig); end
end
