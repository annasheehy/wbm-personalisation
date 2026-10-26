%% fastingFigures — literature-vs-model figures for the fasting chapter
%
% Builds the comparison figures for the fasted Harvey/Harvetta models:
%   FIG 1  Fasted model vs expected 24 h fasting fluxes. Per benchmark: the FASTED
%          L2 (min-norm) flux, the FVA range that flux could have taken under the
%          same conditions, and the literature band. Every value is ONE named
%          reaction's flux (never a sum); the reaction id is printed under each
%          label and the full provenance goes to fig1_litVsModel_provenance.xlsx.
%   FIG 2  Organ x fuel matrix: the fasted exchange flux of every organ for the key
%          fuels, on a signed (uptake/release) diverging scale, with the cells that
%          CHANGED DIRECTION between fed and fasting outlined.
%   FIG 3  Whole-body respiratory quotient, EX_co2[a] / -EX_o2[a] (a ratio, not a flux).
%
% WHY FVA AS WELL AS THE QP POINT. The L2 QP returns the one flux distribution the
% model chose; the nominal lb/ub say nothing about what the rest of the network
% actually permits. fastingTargetFVA minimises and maximises each target reaction
% under the SAME conditions as the BMR solve (brain >= 3500, heart >= 1000,
% objective = 1), so the bar answers "could the model have matched the literature?"
% and the dot answers "did it?".
%
% Units throughout: mmol/person/day. Flux sign convention: positive = release from
% organ into blood/urine, negative = uptake. This is verified per reaction from the
% stoichiometry rather than assumed (see buildFigData).
%
% The models are steady state, so these fluxes are DAY-AVERAGED rates and the
% literature comparator is the average over an overnight-to-24 h fast. See
% fastingLiteratureTargets.m for the bands and their provenance.
%
% Re-runnable on its own: caches ckpt_figData.mat (bounds + fluxes) and
% ckpt_fvaTargets.mat (FVA), so only the first run pays the cost.
%
% Anna Sheehy 2026.

P = persWBMpaths();
resPath = fullfile(P.results, 'fasting');
figPath = fullfile(resPath, 'figures');
if ~exist(figPath, 'dir'), mkdir(figPath); end
figData = fullfile(resPath, 'ckpt_figData.mat');
fvaData = fullfile(resPath, 'ckpt_fvaTargets.mat');

% Okabe-Ito, colourblind-safe: blue = male, vermillion = female
COL = struct('male', [0 114 178]/255, 'female', [213 94 0]/255, ...
             'grey', [115 115 115]/255, 'ink', [0.15 0.15 0.15], 'band', [0.85 0.87 0.90]);
% Three shades per sex for FVA min / QP / FVA max. Within a sex the hue is constant and
% only lightness changes, which stays separable under colour-vision deficiency; the QP
% value is the darkest and largest marker because it is the one the model chose.
SHADE = struct( ...
    'male',   [158 202 225; 8 69 148; 66 146 198]/255, ...      % [min; QP; max]
    'female', [252 174 145; 152 12 19; 222 45 38]/255);
FS = 9;   % point size for a figure placed at \textwidth (~17 cm) with NO scaling
sexes = {'male','female'};


%% ---------------------------------------------------------------- prep / cache
if exist(figData, 'file')
    load(figData, 'D');
    fprintf('figData loaded from cache\n');
else
    D = buildFigData(resPath);
    save(figData, 'D');
end

Tt = fastingLiteratureTargets([], 'male');
if exist(fvaData, 'file')
    load(fvaData, 'fvaM', 'fvaF');
    fprintf('FVA loaded from cache\n');
else
    % needs COBRA + Gurobi; ~1 min per reaction on the 80k-reaction WBM
    fprintf('No FVA cache - computing (this is the slow step).\n');
    changeCobraSolver('gurobi', 'LP', 0);
    mdl = load(fullfile(resPath,'ckpt_feasibleModels.mat'), 'maleFast', 'femaleFast');
    rl  = unique(cellstr(Tt.rxn(strlength(Tt.rxn) > 0)), 'stable');
    fvaM = fastingTargetFVA(mdl.maleFast,   rl);
    fvaF = fastingTargetFVA(mdl.femaleFast, rl);
    fvaConditions = 'fasted models, brain>=3500, heart>=1000, Whole_body_objective_rxn=1';
    save(fvaData, 'fvaM', 'fvaF', 'fvaConditions');
end
FVA = struct('male', fvaM, 'female', fvaF);


%% =============================================================== FIG 1
% Agreement ratio = (fasted flux, oriented so positive = the expected direction)
% divided by the geometric midpoint of the literature band. 1 = agreement, and
% NEGATIVE means the model runs the exchange the opposite way.
useRow = Tt.fluxComparable & ~isnan(Tt.litLow) & strlength(Tt.rxn) > 0;
rows = find(useRow);

A = struct([]);
for s = 1:2
    T = fastingLiteratureTargets([], sexes{s});
    d = D.(sexes{s});  f = FVA.(sexes{s});
    for k = 1:numel(rows)
        r = rows(k);
        i = find(strcmp(d.rxns, char(T.rxn(r))), 1);
        if isempty(i)
            warning('fastingFigures: reaction %s not in model - row skipped', T.rxn(r));
            continue
        end
        sgn = 1; if strcmp(T.expectDir(r), 'uptake'), sgn = -1; end
        if ~d.posIsRelease(i), sgn = -sgn; end          % reaction written the other way round
        expMid = sqrt(T.litLow(r) * T.litHigh(r));
        j = find(f.rxn == T.rxn(r), 1);
        fv = sort(sgn * [f.fvaMin(j) f.fvaMax(j)]) / expMid;

        a = struct();
        a.label  = char(T.label(r));
        if strlength(T.source(r)) == 0, a.label = [a.label ' ' char(8224)]; end
        a.rxn    = char(T.rxn(r));
        a.sex    = sexes{s};
        a.valFed = d.vFed(i);   a.valFast = d.vFast(i);
        a.litLow = T.litLow(r); a.litHigh = T.litHigh(r);
        a.source = char(T.source(r));  a.nativeUnit = char(T.nativeUnit(r));
        a.nativeLow = T.nativeLow(r);  a.nativeHigh = T.nativeHigh(r);
        a.lbFast = d.lbFast(i);  a.ubFast = d.ubFast(i);
        a.fvaMin = f.fvaMin(j);  a.fvaMax = f.fvaMax(j);
        a.sgn = sgn;                       % +1 if a positive flux is the expected direction
        a.ratioFast = sgn * d.vFast(i) / expMid;
        a.bandLo = T.litLow(r)/expMid;  a.bandHi = T.litHigh(r)/expMid;
        a.fvaLo  = fv(1);  a.fvaHi = fv(2);
        % Agreement by construction: if FVA cannot leave the literature band, the model
        % was never free to disagree (the band is imposed as a bound). Detected from the
        % model itself, not asserted in the target list.
        a.imposedAsBound = ~isnan(a.fvaLo) && a.fvaLo >= a.bandLo - 1e-6 && a.fvaHi <= a.bandHi + 1e-6;
        a.verdict = verdictOf(a);
        if isempty(A), A = a; else, A(end+1) = a; end %#ok<SAGROW>
    end
end

% One small panel per benchmark, each on its OWN axis in mmol/person/day (no ratios,
% no log). Per sex, three fasted values as dots: FVA minimum, the L2 (QP) flux the
% model chose, and FVA maximum, each labelled with its number. Shaded = literature.
% Sign is the model's own convention, so a dot on the wrong side of zero is literally
% the exchange running the wrong way. The axis is scaled to the literature band and
% the QP flux; an FVA extreme beyond it is drawn as a triangle at the panel edge,
% still labelled with its true value.
labels = unique({A.label}, 'stable');
nL   = numel(labels);
nCol = 3;  nRow = ceil(nL/nCol);
ptName = {'FVA min', 'QP (L2)', 'FVA max'};

fig1 = figure('Units','centimeters', 'Position',[2 2 17 4.6*nRow+2.8], 'Color','w');
tl1  = tiledlayout(nRow, nCol, 'TileSpacing','compact', 'Padding','compact');
for k = 1:nL
    axk = nexttile(tl1); hold(axk,'on');
    m = A(strcmp({A.label}, labels{k}));

    xs = 0;
    for j = 1:numel(m)
        xs = [xs sort(m(j).sgn * [m(j).litLow m(j).litHigh]) m(j).valFast]; %#ok<AGROW>
    end
    pad = 0.30 * max(max(xs) - min(xs), eps);
    xl  = [min(min(xs) - pad, 0), max(max(xs) + pad, 0)];
    edgeTol = 0.22 * diff(xl);              % labels closer than this are stacked

    for j = 1:numel(m)
        y  = 2 - (j-1);                     % male on top, female below
        sh = SHADE.(m(j).sex);
        bnd = sort(m(j).sgn * [m(j).litLow m(j).litHigh]);
        hBand = patch(axk, [bnd(1) bnd(2) bnd(2) bnd(1)], y+[-0.30 -0.30 0.30 0.30], ...
            COL.band, 'EdgeColor','none');

        vals = [m(j).fvaMin, m(j).valFast, m(j).fvaMax];
        inset = 0.04 * diff(xl);            % off-panel markers sit just inside the edge
        xp = min(max(vals, xl(1) + inset), xl(2) - inset);
        % stack labels of neighbouring dots so no number sits on another
        [~, ord] = sort(xp);  lvl = zeros(1,3);
        for p = 2:3
            if xp(ord(p)) - xp(ord(p-1)) < edgeTol, lvl(ord(p)) = lvl(ord(p-1)) + 1; end
        end
        for p = [1 3 2]                     % QP drawn last, on top
            if isnan(vals(p)), continue; end
            if p == 2                       % the QP (L2) flux the model chose
                mk = 'x';
            elseif vals(p) < xl(1), mk = '<';
            elseif vals(p) > xl(2), mk = '>';
            else,                   mk = 'o';
            end
            if p == 2
                plot(axk, xp(p), y, mk, 'MarkerEdgeColor', sh(p,:), 'MarkerSize', 5.5, 'LineWidth', 1.3);
            else
                plot(axk, xp(p), y, mk, 'MarkerFaceColor', sh(p,:), 'MarkerEdgeColor', 'w', ...
                    'MarkerSize', 5.5, 'LineWidth', 0.6);
            end
            fw = 'normal';  if p == 2, fw = 'bold'; end
            txtCol = sh(p,:);  if p ~= 2, txtCol = sh(p,:) * 0.80; end   % same hue, legible on white
            ha = 'center';                  % keep edge labels inside the panel
            if mk == '<', ha = 'left'; elseif mk == '>', ha = 'right'; end
            text(axk, xp(p), y + 0.16 + 0.19*lvl(p), fluxStr(vals(p)), 'Color', txtCol, ...
                'FontSize', FS-3.5, 'FontWeight', fw, ...
                'HorizontalAlignment', ha, 'VerticalAlignment','bottom');
        end
    end
    plot(axk, [0 0], [0.55 2.95], '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.75);
    set(axk, 'YTick', [1 2], 'YTickLabel', {'F','M'}, 'XLim', xl, 'YLim', [0.55 2.95], ...
        'FontSize', FS-2.5, 'TickDir','out', 'Box','off', 'YColor', COL.ink, 'XColor', COL.ink);
    ttl = labels{k};
    if all([m.imposedAsBound]), ttl = [ttl ' *']; end   % agreement is by construction
    title(axk, {ttl, ['\fontsize{6}\color[rgb]{0.45 0.45 0.45}\it' strrep(m(1).rxn,'_','\_')]}, ...
        'FontSize', FS-1.5, 'FontWeight','bold');
    axk.Toolbar.Visible = 'off';
end
title(tl1, 'Fasted model vs expected fluxes, overnight to 24 h fast', 'FontSize', FS+1, 'FontWeight','bold');
xlabel(tl1, {'mmol/person/day    (positive = release into blood, negative = uptake)', ...
    [char(9664) ' ' char(9654) '  FVA extreme beyond the panel edge (true value printed)'], ...
    '*  literature value is also imposed as a model bound, so agreement is by construction'}, 'FontSize', FS-1.5);

% legend: proxy handles, three shades per sex + the literature band
hL = gobjects(1,7);  lt = cell(1,7);  q = 0;
for s = 1:2
    for p = 1:3
        q = q + 1;
        if p == 2
            hL(q) = plot(axk, NaN, NaN, 'x', 'MarkerEdgeColor', SHADE.(sexes{s})(p,:), ...
                'MarkerSize', 5.5, 'LineWidth', 1.3);
        else
            hL(q) = plot(axk, NaN, NaN, 'o', 'MarkerFaceColor', SHADE.(sexes{s})(p,:), ...
                'MarkerEdgeColor','w', 'MarkerSize', 5.5);
        end
        lt{q} = sprintf('%s %s', upper(sexes{s}(1)), ptName{p});
    end
end
hL(7) = hBand;  lt{7} = 'literature, overnight to 24 h';
lgd = legend(hL([1 4 2 5 3 6 7]), lt([1 4 2 5 3 6 7]), 'NumColumns', 4, 'FontSize', FS-1.5, 'Box','off');
lgd.ItemTokenSize = [10 8];
lgd.Layout.Tile = 'south';

exportgraphics(fig1, fullfile(figPath,'fig1_litVsModel.png'), 'Resolution', 300);
exportgraphics(fig1, fullfile(figPath,'fig1_litVsModel.pdf'), 'ContentType','vector');

% ---- provenance table: exactly which reaction (never a sum) each value comes from
Tp = table(string({A.sex})', string({A.label})', string({A.rxn})', ...
    arrayfun(@(a) sprintf('%g-%g %s', a.nativeLow, a.nativeHigh, a.nativeUnit), A, 'UniformOutput', false)', ...
    [A.litLow]', [A.litHigh]', string({A.source})', [A.valFed]', [A.valFast]', ...
    [A.lbFast]', [A.ubFast]', [A.fvaMin]', [A.fvaMax]', [A.ratioFast]', [A.imposedAsBound]', ...
    string({A.verdict})', ...
    'VariableNames', {'sex','benchmark','modelReaction','literatureAsReported','litLow_mmolPerDay', ...
    'litHigh_mmolPerDay','source','fedFlux','fastedFlux','fastedLb','fastedUb','fvaMin','fvaMax', ...
    'fastedRatio','imposedAsBound','verdict'});
iR = find(strcmp(Tt.id, 'rq'), 1);
for s = 1:2          % RQ: the one derived value = ratio of two single reactions
    d = D.(sexes{s});
    Tp = [Tp; {string(sexes{s}), "Whole-body respiratory quotient", "EX_co2[a] / -EX_o2[a]", ...
        sprintf('%g-%g ratio', Tt.nativeLow(iR), Tt.nativeHigh(iR)), NaN, NaN, string(Tt.source(iR)), ...
        d.vFed(d.iCO2)/-d.vFed(d.iO2), d.vFast(d.iCO2)/-d.vFast(d.iO2), NaN, NaN, NaN, NaN, NaN, false, ...
        sprintf('O2 %.0f (bounds %.0f to %.0f), CO2 %.0f (bounds %.0f to %.0f)', d.vFast(d.iO2), ...
        d.lbFast(d.iO2), d.ubFast(d.iO2), d.vFast(d.iCO2), d.lbFast(d.iCO2), d.ubFast(d.iCO2))}]; %#ok<AGROW>
end
try
    writetable(Tp, fullfile(figPath, 'fig1_litVsModel_provenance.xlsx'), 'WriteMode', 'replacefile');
catch err
    % typically the workbook is open in Excel - do not lose the figures over it
    warning('fastingFigures:xlsxLocked', 'Could not write the provenance xlsx (%s). Close it and re-run.', err.message);
end
disp(Tp(:, {'sex','benchmark','modelReaction','fastedFlux','fvaMin','fvaMax','verdict'}));


%% =============================================================== FIG 2
% Organ x fuel matrix of the FASTED flux, signed (release positive), on an asinh
% (symmetric log) colour scale so that both ~1 and ~1000 mmol/day are readable.
fuels = {'glc_D','lac_L','pyr','glyc','bhb','acac','ala_L','gln_L','hdca','ocdca','ocdcea','urea','o2','co2'};
fuelLab = {'glucose','lactate','pyruvate','glycerol','3-OHB','AcAc','alanine','glutamine', ...
           'palmitate','stearate','oleate','urea','O_2','CO_2'};
organs = {'Liver','Adipocytes','Muscle','Brain','Heart','Kidney','RBC','sIEC','Colon','Lung','Skin'};

fig2 = figure('Units','centimeters', 'Position',[2 2 17 12], 'Color','w');
tl2 = tiledlayout(2, 1, 'TileSpacing','compact', 'Padding','compact');
for s = 1:2
    d = D.(sexes{s});
    Mfast = nan(numel(organs), numel(fuels));  flip = false(size(Mfast));
    for o = 1:numel(organs)
        for f = 1:numel(fuels)
            i = findExch(d, organs{o}, fuels{f});
            if isempty(i), continue; end
            sg = 1; if ~d.posIsRelease(i), sg = -1; end
            Mfast(o,f) = sg * d.vFast(i);
            % a reversal needs BOTH fluxes clearly non-zero (not "stopped" / solver noise)
            flip(o,f)  = sign(d.vFed(i)) * sign(d.vFast(i)) < 0 && min(abs([d.vFed(i) d.vFast(i)])) > 0.5;
        end
    end
    ax = nexttile(tl2, s);
    imagesc(ax, asinh(Mfast/2), 'AlphaData', ~isnan(Mfast));  hold(ax,'on');
    set(ax, 'Color', [0.80 0.80 0.80]);  axis(ax, 'tight');     % grey = organ has no such exchange
    ax.Toolbar.Visible = 'off';
    caxis(ax, asinh([-1 1]*3000/2));
    colormap(ax, divergingMap(256));
    for o = 1:numel(organs)
        for f = 1:numel(fuels)
            if flip(o,f)
                rectangle(ax, 'Position', [f-0.5 o-0.5 1 1], 'EdgeColor', COL.ink, 'LineWidth', 1.4);
            end
        end
    end
    set(ax, 'XTick', 1:numel(fuels), 'XTickLabel', fuelLab, 'YTick', 1:numel(organs), ...
        'YTickLabel', organs, 'FontSize', FS-1, 'TickDir','out', 'TickLength',[0 0]);
    xtickangle(ax, 45);
    names = {'Harvey (male)', 'Harvetta (female)'};
    title(ax, sprintf('%s   %s, fasted', char(64+s), names{s}), 'FontSize', FS+1, 'FontWeight','bold');
    if s == 2
        cb = colorbar(ax, 'southoutside');
        tks = [-3000 -300 -30 0 30 300 3000];
        set(cb, 'Ticks', asinh(tks/2), 'TickLabels', arrayfun(@(x) sprintf('%d',x), tks, 'UniformOutput',false), ...
            'FontSize', FS-1.5);
        cb.Label.String = 'mmol/person/day   red = release into blood, blue = uptake;  outlined = reversed vs fed;  grey = no exchange';
        cb.Label.FontSize = FS-1;
    end
end
exportgraphics(fig2, fullfile(figPath,'fig2_organFuelMatrix.png'), 'Resolution', 300);
exportgraphics(fig2, fullfile(figPath,'fig2_organFuelMatrix.pdf'), 'ContentType','vector');


%% =============================================================== FIG 3
% Respiratory quotient: a ratio, not a flux, so it gets its own small figure.
fig3 = figure('Units','centimeters', 'Position',[2 2 7 7], 'Color','w');
ax2 = axes(fig3); hold(ax2,'on');
patch(ax2, [0.6 2.4 2.4 0.6], [Tt.nativeLow(iR) Tt.nativeLow(iR) Tt.nativeHigh(iR) Tt.nativeHigh(iR)], ...
    COL.band, 'EdgeColor','none');
yline(ax2, 0.70, ':', 'Color', COL.ink, 'LineWidth', 1, 'FontSize', FS-1, ...
    'Label','0.70 pure fat', 'LabelHorizontalAlignment','right', 'LabelVerticalAlignment','top');
for s = 1:2
    d = D.(sexes{s});
    rqFed  = d.vFed(d.iCO2)  / -d.vFed(d.iO2);
    rqFast = d.vFast(d.iCO2) / -d.vFast(d.iO2);
    plot(ax2, [1 2], [rqFed rqFast], '-', 'Color', COL.(sexes{s}), 'LineWidth', 1.2);
    plot(ax2, 1, rqFed,  'o', 'MarkerEdgeColor', COL.(sexes{s}), 'MarkerFaceColor','w', 'MarkerSize',6, 'LineWidth',1.2);
    plot(ax2, 2, rqFast, 'o', 'MarkerEdgeColor', COL.(sexes{s}), 'MarkerFaceColor', COL.(sexes{s}), 'MarkerSize',6, 'LineWidth',1.2);
    fprintf('RQ %-6s fed %.3f -> fasted %.3f\n', sexes{s}, rqFed, rqFast);
end
set(ax2, 'XTick', [1 2], 'XTickLabel', {'fed','fasted'}, 'FontSize', FS, 'TickDir','out', ...
    'Box','off', 'YColor', COL.ink, 'XColor', COL.ink);
xlim(ax2, [0.6 2.4]); ylim(ax2, [0.4 1.0]);
ylabel(ax2, 'respiratory quotient', 'FontSize', FS);
title(ax2, {'Whole-body RQ', '\rm\fontsize{7}EX\_co2[a] / -EX\_o2[a]'}, 'FontSize', FS+1, 'FontWeight','bold');
ax2.Toolbar.Visible = 'off';
exportgraphics(fig3, fullfile(figPath,'fig3_RQ.png'), 'Resolution', 300);
exportgraphics(fig3, fullfile(figPath,'fig3_RQ.pdf'), 'ContentType','vector');

fprintf('\nfigures written to %s\n', figPath);


%% =============================================================== helpers
function D = buildFigData(resPath)
% Pull the reaction names, bounds and L2 fluxes out of the big checkpoints once.
% posIsRelease records, PER REACTION and from the stoichiometry, whether a positive
% flux means release out of the organ (into [bc]/[csf]/[u]) - never assumed.
t0 = tic;
load(fullfile(resPath,'ckpt_feasibleModels.mat'), 'maleFed','femaleFed','maleFast','femaleFast');
load(fullfile(resPath,'ckpt_fluxL2.mat'), 'fluxFedM','fluxFedF','fluxFastM','fluxFastF');
sets = {'male', maleFed, maleFast, fluxFedM, fluxFastM; 'female', femaleFed, femaleFast, fluxFedF, fluxFastF};
for s = 1:2
    mFed = sets{s,2};  mFast = sets{s,3};
    d = struct();
    d.rxns  = mFast.rxns;
    d.lbFed = mFed.lb;   d.ubFed = mFed.ub;
    d.lbFast = mFast.lb; d.ubFast = mFast.ub;
    d.vFed  = full(sets{s,4});  d.vFast = full(sets{s,5});
    d.iO2   = find(strcmp(d.rxns,'EX_o2[a]'), 1);
    d.iCO2  = find(strcmp(d.rxns,'EX_co2[a]'), 1);
    % orientation: coefficient of the OUTSIDE metabolite ([bc]/[csf]/[u]/[a])
    n = numel(d.rxns);
    d.posIsRelease = true(n,1);
    isEx = ~cellfun(@isempty, regexp(d.rxns, '_EX_|^EX_', 'once'));
    for i = find(isEx(:)')
        mr  = find(mFast.S(:,i) ~= 0);
        out = mr(~cellfun(@isempty, regexp(mFast.mets(mr), '\[(bc|csf|u|a)\]$', 'once')));
        if isempty(out), continue; end
        d.posIsRelease(i) = full(mFast.S(out(1), i)) > 0 || ...
                            (numel(mr) == 1 && full(mFast.S(mr,i)) < 0);   % sink: met[u] <=>
    end
    nRev = sum(~d.posIsRelease);
    if nRev > 0
        fprintf('%s: %d exchange reaction(s) written with uptake positive - orientation corrected per reaction.\n', ...
            sets{s,1}, nRev);
    end
    D.(sets{s,1}) = d;
end
fprintf('buildFigData: %.0f s\n', toc(t0));
end

function v = verdictOf(a)
% Plain-words reading of one row: how the fasted flux compares with the literature,
% and whether the literature is reachable at all (from FVA, not from lb/ub).
if a.ratioFast < -1e-3
    match = 'reversed';
elseif abs(a.ratioFast) <= 1e-3
    match = 'essentially zero';
elseif a.ratioFast >= a.bandLo*(1-1e-6) && a.ratioFast <= a.bandHi*(1+1e-6)
    match = 'agrees';
else
    match = sprintf('%.2g%s expected', a.ratioFast, char(215));
end
if isnan(a.fvaLo)
    reach = 'FVA failed';
elseif a.fvaHi - a.fvaLo < 1e-6*max(1, abs(a.fvaHi))
    reach = 'flux is fixed';
elseif a.fvaHi >= a.bandLo && a.fvaLo <= a.bandHi
    reach = 'reachable';
else
    reach = 'UNREACHABLE';
end
v = [match ', ' reach];
end

function s = fluxStr(x)
if abs(x) < 0.01, s = sprintf('%.1g', x); else, s = sprintf('%.0f', x); end
end

function i = findExch(d, organ, met)
i = find(strcmp(d.rxns, sprintf('%s_EX_%s(e)_[bc]', organ, met)), 1);
if isempty(i)
    i = find(strcmp(d.rxns, sprintf('%s_EX_%s(e)_[csf]', organ, met)), 1);
end
end

function cm = divergingMap(n)
% Two hues with a neutral midpoint (never a rainbow): blue = uptake, red = release.
lo = [ 33  102 172]/255;  mid = [0.96 0.96 0.96];  hi = [178  24  43]/255;
h = floor(n/2);
cm = [ [linspace(lo(1),mid(1),h)' linspace(lo(2),mid(2),h)' linspace(lo(3),mid(3),h)'] ; ...
       [linspace(mid(1),hi(1),n-h)' linspace(mid(2),hi(2),n-h)' linspace(mid(3),hi(3),n-h)'] ];
end
