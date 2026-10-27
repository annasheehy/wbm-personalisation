function pct = metabFeasibilityFigure(feas, fixTest, windows, figPath, thesisFig)
% Figure of metabolomicsPersonalisation.m (metabC_feasibility): feasible personalised models
% (%) at each measured concentration allowance in Harvey (left) and Harvetta (right), as
% built by persWBMmetabolomics and with the measurement given priority where a personalised
% bound crossed a fixed bound. Feasibility with the resting brain and heart ATP demands
% (feasibleFloors) was identical and is returned but not drawn.
%
% INPUT:  feas     one row per sex, profile and allowance (feasibility.csv)
%         fixTest  models with crossed bounds, retested with the measurement given priority (fixTest.csv)
% OUTPUT: pct.(sex) = nWindows x 4 [as built, + ATP demands, priority, priority + ATP]
%
% Anna Sheehy 2026.

sexes = {'male', 'female'};  sexLab = {'Harvey', 'Harvetta'};
colours = [0 114 178; 0 158 115] / 255;                         % Okabe-Ito blue, green
serLab = {'As built', 'Measurement given priority'};
labW = compose('\\pm%d%%', round(100 * windows));
pct = struct();
h = figure('Units', 'centimeters', 'Position', [2 2 16 6.5], 'Color', 'w');
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for k = 1:2
    s = sexes{k};  P = zeros(numel(windows), 4);
    for wI = 1:numel(windows)
        r  = feas(strcmp(feas.sex, s) & abs(feas.window - windows(wI)) < 1e-9, :);
        fx = fixTest(strcmp(fixTest.sex, s) & abs(fixTest.window - windows(wI)) < 1e-9, :);
        P(wI, :) = 100 * [mean(r.feasible), mean(r.feasibleFloors), ...
            mean(r.feasible | ismember(r.individual, fx.individual(fx.feasible == 1))), ...
            mean(r.feasibleFloors | ismember(r.individual, fx.individual(fx.feasibleFloors == 1)))];
    end
    pct.(s) = P;
    nexttile;
    b = bar(categorical(labW, labW), P(:, [1 3]), 'BarWidth', 0.7);
    for j = 1:2, b(j).FaceColor = colours(j, :); end
    set(gca, 'FontSize', 8, 'TickDir', 'out', 'Box', 'off', 'YLim', [0 105]);
    if k == 1, ylabel('Feasible models (%)'); else, set(gca, 'YTickLabel', []); end
    xlabel('Measured concentration allowance');
    title(sexLab{k}, 'FontWeight', 'normal');
    if k == 2
        lg = legend(b, serLab, 'NumColumns', 2, 'Box', 'off');
        lg.Layout.Tile = 'south';
    end
end
exportgraphics(h, fullfile(figPath, 'metabC_feasibility.png'), 'Resolution', 300);
exportgraphics(h, fullfile(figPath, 'metabC_feasibility.pdf'), 'ContentType', 'vector');
if exist(thesisFig, 'dir'), copyfile(fullfile(figPath, 'metabC_feasibility.png'), thesisFig); end
end
