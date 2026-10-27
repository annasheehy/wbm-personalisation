%% supplementaryMaterial.m
% Thesis Ch5: supplementary material from the saved inputs and results (no solves).
%   ThesisChapter/LateX/appendix_chapter5.tex   Tables S1-S9 for the thesis appendix
%   results/supplementary/SupplementaryData_Ch5_Inputs.xlsx    healthy values, ranges, benchmarks, cohort
%   results/supplementary/SupplementaryData_Ch5_Results.xlsx   Tables S4-S9 and the full results
% Run after the other scripts (README run order); needs oatBMRdecomposition.csv and
% fastingSetupTests.csv.
%
% Anna Sheehy 2026.

P = persWBMpaths();
fastPath = fullfile(P.results, 'fasting');  metPath = fullfile(P.results, 'metabolomics');
outPath  = fullfile(P.results, 'supplementary');  if ~exist(outPath, 'dir'), mkdir(outPath); end
sexes = {'male', 'female'};  sexLab = {'Harvey', 'Harvetta'};
IPd  = load(fullfile(fastPath, 'ckpt_HMDBconstrained.mat'), 'ipMale', 'ipFemale');
ipar = struct('male', IPd.ipMale, 'female', IPd.ipFemale);
X = struct();                                             % tables for the workbooks
tex = {};                                                 % LaTeX of the appendix

%% S1 healthy blood metabolite values -----------------------------------------------------
refV  = readtable(fullfile(P.references, 'healthyReferenceValues.csv'), 'TextType', 'char', 'Delimiter', ',');
panel = healthyReferencePanel(fullfile(P.references, 'healthyReferenceValues.csv'));
inPanel = readtable(fullfile(metPath, 'inputPanel.csv'), 'TextType', 'char', 'Delimiter', ',');
refV = refV(ismember(refV.vmh, inPanel.vmh), :);
rows = cell(height(refV), 6);
for r = 1:height(refV)
    q = strcmp(panel.vmh, refV.vmh{r});
    switch refV.sex{r}
        case 'female', fit = [panel.p025Female(q) panel.p975Female(q)];  sx = 'F';
        case 'male',   fit = [panel.p025Male(q) panel.p975Male(q)];      sx = 'M';
        otherwise,     fit = [panel.p025Male(q) panel.p975Male(q)];      sx = 'M, F';
    end
    [st, rep] = reportedStat(refV.stat{r}, refV.v1(r), refV.v2(r), refV.v3(r));
    rows(r, :) = {refV.name{r}, sx, st, sprintf('%s %s', rep, unitTex(refV.unit{r})), ...
        sprintf('%s--%s', fmt(fit(1)), fmt(fit(2))), sprintf('\\cite{%s}', refV.source{r})};
end
tex{end+1} = latexTable('Healthy blood concentrations as reported and the fitted distributions', ...
    ['Values as reported by each source and the 2.5th to 97.5th percentile of the log-normal distribution fitted to them, from which the synthetic profiles were drawn (Methods). ' ...
     'M, F: one distribution for both sexes, because the source reports men and women together.'], ...
    'tab:suppBlood', {'Metabolite', 'Sex', 'Statistic', 'Reported', 'Fitted ($\mu$M)', 'Source'}, 'llllrl', rows);
X.in.S1_bloodValues = [refV, table(arrayfun(@(r) S1fit(panel, refV, r, 1), (1:height(refV))'), ...
    arrayfun(@(r) S1fit(panel, refV, r, 2), (1:height(refV))'), 'VariableNames', {'fitted_p025_uM', 'fitted_p975_uM'})];
X.in.S1_fittedPanel = panel;

%% S2 physiological ranges ----------------------------------------------------------------
refP = readtable(fullfile(P.references, 'healthyPhysiologyRanges.csv'), 'TextType', 'char', 'Delimiter', ',');
[~, physTab] = healthyPhysiologyRanges(fullfile(P.references, 'healthyPhysiologyRanges.csv'), ipar);
inLab = containers.Map({'bodyWeight', 'Height', 'HeartRate', 'Hematocrit', 'GlomerularFiltrationRate', 'UrineCreatinine'}, ...
    {'Body mass index', 'Height', 'Resting heart rate', 'Haematocrit', 'GFR', 'Urine creatinine'});
order = {'bodyWeight', 'Height', 'HeartRate', 'Hematocrit', 'GlomerularFiltrationRate', 'UrineCreatinine'};
[~, o] = sort(cellfun(@(f) find(strcmp(order, f)), refP.field));  refP = refP(o, :);
rows = cell(height(refP), 6);
for r = 1:height(refP)
    age = regexp(refP.note{r}, 'age (\d+)-(\d+) y', 'tokens', 'once');
    if isempty(age), grp = '--'; else, grp = sprintf('%s--%s', age{1}, age{2}); end
    [st, rep] = reportedStat(refP.stat{r}, refP.v1(r), refP.v2(r), refP.v3(r));
    if strcmp(refP.field{r}, 'bodyWeight'), st = 'normal range'; end
    nn = '';  if ~isnan(refP.n(r)), nn = sprintf(' (n = %d)', refP.n(r)); end
    rows(r, :) = {inLab(refP.field{r}), upper(refP.sex{r}(1)), grp, st, ...
        [strtrim(sprintf('%s %s', rep, unitTex(refP.unit{r}))) nn], sprintf('\\cite{%s}', refP.source{r})};
end
tex{end+1} = latexTable('Physiological ranges as reported', ...
    ['Values as reported by each source, by sex (M, F) and age group (years). Rows of one input and sex were pooled by sample size, ' ...
     'and the ranges tested were calculated from them and converted to the units of the models (Methods, \autoref{tab:physRanges}).'], ...
    'tab:suppPhys', {'Input', 'Sex', 'Age', 'Statistic', 'Reported', 'Source'}, 'llllll', rows);
X.in.S2_physiologyValues = refP;  X.in.S2_rangesTested = physTab;

%% S3 fasting benchmarks ------------------------------------------------------------------
Tm = fastingLiteratureTargets([], 'male');  Tf = fastingLiteratureTargets([], 'female');
use = find(Tm.fluxComparable & ~isnan(Tm.litLow));
rows = cell(numel(use), 6);
for k = 1:numel(use)
    r = use(k);  lab = char(Tm.label(r));
    if ~Tm.verified(r), lab = [lab '$^{a}$']; end
    rows(k, :) = {lab, sprintf('%g--%g %s', Tm.nativeLow(r), Tm.nativeHigh(r), unitTex(char(Tm.nativeUnit(r)))), ...
        regexprep(strrep(strrep(char(Tm.window(r)), '->', 'to'), ' (obligate)', ''), '(\d)-(\d)', '$1--$2'), sprintf('%s--%s', fmt(Tm.litLow(r)), fmt(Tm.litHigh(r))), ...
        sprintf('%s--%s', fmt(Tf.litLow(r)), fmt(Tf.litHigh(r))), sprintf('\\cite{%s}', char(Tm.source(r)))};
end
tex{end+1} = latexTable('Fluxes measured in fasting humans, as reported and in model units', ...
    ['Ranges as reported and converted to mmol per day for the body weight of each model (70~kg Harvey, 58~kg Harvetta): ' ...
     '$\mu$mol\,kg$^{-1}$\,min$^{-1}$ $\times$ body weight $\times$ 1.44, and mg\,min$^{-1}$ per 65~kg $\times$ 1.44 $\times$ body weight / 65 / molar mass. ' ...
     '$^{a}$A single reported rate $\pm 20\%$ (red blood cell glucose uptake, muscle and hepatic alanine exchange), or with a lower end 20\% below it (brain glucose uptake).'], ...
    'tab:suppBench', {'Flux', 'Reported', 'Fasting time', 'Harvey', 'Harvetta', 'Source'}, 'lllrrl', rows);
X.in.S3_fastingBenchmarks = [table(repmat("Harvey", height(Tm), 1), 'VariableNames', {'model'}), Tm; ...
                             table(repmat("Harvetta", height(Tf), 1), 'VariableNames', {'model'}), Tf];

%% S4 crossed-bound thresholds (alanine) --------------------------------------------------
thr = readtable(fullfile(metPath, 'conflictThresholds.csv'), 'TextType', 'char', 'Delimiter', ',');
rows = cell(height(thr), 4);
for r = 1:height(thr)
    rows(r, :) = {sexLab{strcmp(sexes, thr.sex{r})}, sprintf('$\\pm %d\\%%$', round(100 * thr.window(r))), ...
        sprintf('%.0f', thr.threshold_uM(r)), sprintf('%.0f', 100 * thr.popFraction(r))};
end
tex{end+1} = latexTable('Blood alanine concentrations at which the hepatic alanine bounds cross', ...
    ['For each measured concentration allowance, the blood alanine concentration below which the personalised limit on hepatic alanine uptake falls short of the fixed requirement of the control model, ' ...
     'and the share of the fitted healthy distribution below it, which is the same distribution for both sexes \cite{Trabado2017TheAge}.'], ...
    'tab:suppThresh', {'Model', 'Allowance', 'Threshold ($\mu$M)', 'Healthy values below (\%)'}, 'lrrr', rows);
X.res.S4_alanineThresholds = thr;

%% S5 BMR change with heart rate and haematocrit, by organ ------------------------------
dec = readtable(fullfile(metPath, 'oatBMRdecomposition.csv'), 'TextType', 'char', 'Delimiter', ',');
cols = {'Heart rate', 'male'; 'Heart rate', 'female'; 'Haematocrit', 'male'; 'Haematocrit', 'female'};
org = dec(isnan(dec.low), :);
organs = unique(org.item);  V = zeros(numel(organs), 4);
for c = 1:4
    sel = strcmp(org.input, cols{c, 1}) & strcmp(org.sex, cols{c, 2});
    [~, ii] = ismember(org.item(sel), organs);  V(ii, c) = org.change(sel);
end
[~, o] = sort(max(abs(V), [], 2), 'descend');  top = o(1:8);
rows = cell(0, 5);
for i = top', rows(end+1, :) = [{organName(organs{i})}, arrayfun(@signed, V(i, :), 'UniformOutput', false)]; end %#ok<SAGROW>
rows(end+1, :) = [{'Other organs'}, arrayfun(@signed, sum(V, 1) - sum(V(top, :), 1), 'UniformOutput', false)];
item = {'BMR (kcal/day)', 'Total (BMR change)', 1; 'EX_o2[a]', 'Whole-body O$_2$ uptake (mmol/day)', -1};
for k = 1:size(item, 1)
    v = nan(1, 4);
    for c = 1:4
        sel = strcmp(dec.item, item{k, 1}) & strcmp(dec.input, cols{c, 1}) & strcmp(dec.sex, cols{c, 2});
        if any(sel), v(c) = item{k, 3} * dec.change(sel); end
    end
    rows(end+1, :) = [item(k, 2), arrayfun(@signed, v, 'UniformOutput', false)]; %#ok<SAGROW>
end
tex{end+1} = latexTable('Change in BMR with heart rate and haematocrit, by organ', ...
    ['Change in the ATP consumption counted in the BMR (kcal per day) from the low to the high value of each input (QP-FBA), by organ, and the change in whole-body O$_2$ uptake. CO$_2$ output stayed at its lower limit in all four cases. HR, heart rate; Hct, haematocrit.'], ...
    'tab:suppBMRorgan', {'', 'HR, Harvey', 'HR, Harvetta', 'Hct, Harvey', 'Hct, Harvetta'}, 'lrrrr', rows);
X.res.S5_bmrByOrgan = dec;

%% S6 height through cardiac output and blood volume -----------------------------------
warning('off', 'MATLAB:load:variableNameIsKeyword');    % saved field 'case' loads as case_
caseOrder = {'std_default', 'pers_default', 'defaultHeight', 'lowHeight', 'highHeight', ...
             'opt4pers_default', 'opt4pers_lowH', 'opt4pers_highH', 'pers_lowHct', 'pers_highHct'};
coLab  = containers.Map({'default', 'opt4'}, {'HR $\times$ SV', 'BSA'});
orgLab = containers.Map({'none', 'std', 'pers'}, {'model', 'body weight', 'personalised BV'});
rows = cell(0, 8);  H = table();
for k = 1:2
    for c = 1:numel(caseOrder)
        f = fullfile(metPath, 'heightCO', sprintf('%s_%s.mat', sexes{k}, caseOrder{c}));
        if ~exist(f, 'file'), continue, end
        R = load(f);
        rows(end+1, :) = {sexLab{k}, sprintf('%.1f', R.height), sprintf('%.3f', R.hct), coLab(R.coMode), orgLab(R.orgMode), ...
            sprintf('%.0f', R.CO), sprintf('%.0f', R.kcal), sprintf('%.0f', -R.o2)}; %#ok<SAGROW>
        H = [H; table(string(sexLab{k}), string(caseOrder{c}), R.height, R.hct, string(R.coMode), string(R.orgMode), R.CO, R.BV, R.kcal, -R.o2, ...
            'VariableNames', {'model', 'case', 'height_cm', 'haematocrit', 'cardiacOutput', 'organWeights', 'CO_mLmin', 'bloodVolume_mL', 'BMR_kcal', 'O2uptake_mmol'})]; %#ok<AGROW>
    end
end
tex{end+1} = latexTable('Height through cardiac output and blood volume', ...
    ['Cardiac output (CO) from heart rate and stroke volume (HR $\times$ SV, the default) or from body surface area (BSA; \texttt{individualizedLabReport}, option 4, in mL per minute). ' ...
     'Organ weights: those of the model, recalculated from body weight (\texttt{calcOrganFract}), or with the personalised blood volume (BV; \texttt{calcOrganFractPersonal}, red blood cell weight scaled by haematocrit). ' ...
     'BMR and O$_2$ uptake from QP-FBA. The control models have a BMR of 1673 (Harvey) and 1426~kcal per day (Harvetta).'], ...
    'tab:suppHeight', {'Model', 'Height (cm)', 'Hct', 'CO from', 'Organ weights', 'CO (mL/min)', 'BMR (kcal/day)', 'O$_2$ uptake (mmol/day)'}, 'lrrllrrr', rows);
X.res.S6_heightCases = H;

%% S7 fasted-state setup tests -----------------------------------------------------------
st = readtable(fullfile(fastPath, 'fastingSetupTests.csv'), 'TextType', 'char', 'Delimiter', ',');
[~, keep] = unique(strcat(st.sex, '_', string(st.aaRelease), '_', string(st.phosphateLimit)), 'stable');
st = st(sort(keep), :);
rows = cell(height(st), 5);
for r = 1:height(st)
    pl = 'open';  if isfinite(st.phosphateLimit(r)), pl = sprintf('%.1f', st.phosphateLimit(r)); end
    rows(r, :) = {sexLab{strcmp(sexes, st.sex{r})}, sprintf('%d', st.aaRelease(r)), pl, ...
        ifNaN(st.brainATPmax(r), 'infeasible'), ifNaN(st.heartATPmax(r), 'infeasible')};
end
tex{end+1} = latexTable('Tests that set the fasted state', ...
    ['Largest possible brain and heart ATP hydrolysis (LP, mmol per day) of the fasted models with whole-body maintenance fixed to one, for each release limit of the amino-acid stores and with phosphate uptake capped at its dietary provision (\texttt{EUAverageDietNew}) or left open. ' ...
     'The resting demands used for BMR are 3500 (brain) and 1000~mmol per day (heart).'], ...
    'tab:suppSetup', {'Model', 'Amino-acid release (mmol/day)', 'Phosphate uptake limit (mmol/day)', 'Brain ATP', 'Heart ATP'}, 'lrrrr', rows);
X.res.S7_fastedSetupTests = st;

%% S8 fed and fasted fluxes of the 11 read-outs -----------------------------------------
Dd = load(fullfile(fastPath, 'ckpt_figData.mat'), 'D');  D = Dd.D;
ro = oatReadouts();
rows = cell(size(ro, 1), 5);  FF = table();
for i = 1:size(ro, 1)
    v = nan(1, 4);
    for k = 1:2
        d = D.(sexes{k});  j = find(strcmp(d.rxns, ro{i, 1}), 1);
        sg = 1;  if ~d.posIsRelease(j), sg = -1; end
        v(2 * k - 1 : 2 * k) = sg * [d.vFed(j) d.vFast(j)];
    end
    rows(i, :) = [{ro{i, 2}}, arrayfun(@fmt, v, 'UniformOutput', false)];
    FF = [FF; table(string(ro{i, 2}), string(ro{i, 1}), v(1), v(2), v(3), v(4), 'VariableNames', ...
        {'exchange', 'reaction', 'HarveyFed', 'HarveyFasted', 'HarvettaFed', 'HarvettaFasted'})]; %#ok<AGROW>
end
tex{end+1} = latexTable('Fed and fasted fluxes of the 11 organ exchanges', ...
    ['Basal flux distribution (QP-FBA), mmol per day; positive values are release from the organ (into the circulating blood, cerebrospinal fluid or urine) and negative values uptake. ' ...
     'Liver and kidney exchanges are those with the circulating blood; kidney glucose is glomerular filtration.'], ...
    'tab:suppFedFast', {'Exchange', 'Harvey fed', 'Harvey fasted', 'Harvetta fed', 'Harvetta fasted'}, 'lrrrr', rows);
X.res.S8_fedFastedFluxes = FF;

%% S9 urine creatinine in the fasted models ---------------------------------------------
M  = load(fullfile(fastPath, 'ckpt_feasibleModels.mat'), 'maleFast', 'femaleFast');
L2 = load(fullfile(fastPath, 'ckpt_fluxL2.mat'), 'fluxFastM', 'fluxFastF');
Fv = load(fullfile(fastPath, 'ckpt_fvaTargets.mat'), 'fvaM', 'fvaF');
mods = {M.maleFast, M.femaleFast};  vv = {L2.fluxFastM, L2.fluxFastF};  fv = {Fv.fvaM, Fv.fvaF};
rows = cell(0, 5);  U = table();
for k = 1:2
    m = mods{k};  iu = strcmp(m.rxns, 'EX_urea[u]');  j = find(fv{k}.rxn == "EX_urea[u]", 1);
    r0 = {sprintf('%.1f--%.1f (default)', ipar.(sexes{k}).MConUrCreatinineMin, ipar.(sexes{k}).MConUrCreatinineMax), 'yes', ...
          bmrFromFlux(m, vv{k}), full(vv{k}(iu)), [fv{k}.fvaMin(j) fv{k}.fvaMax(j)]};
    allR = {r0};
    for L = 1:2
        Rr = load(fullfile(fastPath, 'urineCreatinine', sprintf('%s_%d.mat', sexes{k}, L)), 'r');  r = Rr.r;
        allR{end+1} = {sprintf('%.1f', r.uc), ifelseStr(r.feasStat == 1, 'yes', 'no'), r.kcal, r.ureaQP, r.ureaFVA}; %#ok<SAGROW>
    end
    for a = 1:numel(allR)
        x = allR{a};
        rows(end+1, :) = {sexLab{k}, x{1}, x{2}, sprintf('%.0f', x{3}), sprintf('%.2f (%.2f to %s)', x{4}, x{5}(1), fmt(x{5}(2)))}; %#ok<SAGROW>
        U = [U; table(string(sexLab{k}), string(x{1}), string(x{2}), x{3}, x{4}, x{5}(1), x{5}(2), 'VariableNames', ...
            {'model', 'urineCreatinine_mgdL', 'feasible', 'BMR_kcal', 'ureaQP', 'ureaFVAmin', 'ureaFVAmax'})]; %#ok<AGROW>
    end
end
tex{end+1} = latexTable('Urine creatinine in the fasted models', ...
    ['Fasted models with the default urine creatinine concentration and with the low and high end of the range tested (\autoref{tab:physRanges}). ' ...
     'BMR and urinary urea excretion from QP-FBA, with the feasible range of urea excretion (FVA) in brackets, in mmol per day; measured excretion after one day of fasting is 180 to 220~mmol per day \cite{Giesecke1989ProteinMethylhistidines}.'], ...
    'tab:suppUrCr', {'Model', 'Urine creatinine (mg/dL)', 'Feasible', 'BMR (kcal/day)', 'Urinary urea (mmol/day)'}, 'llrrl', rows);
X.res.S9_urineCreatinineFasted = U;

%% full results for the workbook ---------------------------------------------------------
rd = @(f) readtable(f, 'TextType', 'char', 'Delimiter', ',');
X.res.oatEffects      = rd(fullfile(metPath, 'oatEffects_long.csv'));
X.res.oatRanking      = rd(fullfile(metPath, 'oatRanking.csv'));
X.res.oatL2Ranking    = rd(fullfile(metPath, 'oatL2Ranking.csv'));
X.res.feasibility     = rd(fullfile(metPath, 'feasibility.csv'));
X.res.crossedBounds   = rd(fullfile(metPath, 'conflicts.csv'));
X.res.measurementPriority = rd(fullfile(metPath, 'fixTest.csv'));
X.res.uptakeLimitRatio    = rd(fullfile(metPath, 'boundRatio.csv'));
X.res.activeState     = rd(fullfile(P.results, 'active', 'activeSummary.csv'));
Mp = load(fullfile(metPath, 'ckpt_map.mat'), 'map');
X.res.boundMap_Harvey   = Mp.map.male;  X.res.boundMap_Harvetta = Mp.map.female;
C = load(fullfile(metPath, 'ckpt_cohort.mat'), 'conc', 'metIDs');
ids = [compose('M%02d', 1:size(C.conc.male, 1)), compose('F%02d', 1:size(C.conc.female, 1))]';
X.in.syntheticCohort_uM = [table(string(ids), [repmat("male", size(C.conc.male, 1), 1); repmat("female", size(C.conc.female, 1), 1)], ...
    'VariableNames', {'ID', 'sex'}), array2table([C.conc.male; C.conc.female], 'VariableNames', C.metIDs)];

%% write ----------------------------------------------------------------------------------
desc = struct( ...
  'S1_bloodValues', 'Table S1: healthy blood concentrations as reported, with the fitted 2.5th-97.5th percentiles (uM)', ...
  'S1_fittedPanel', 'log-normal parameters (log uM) and percentiles of the 26 panel metabolites, per sex', ...
  'S2_physiologyValues', 'Table S2: physiological values as reported', ...
  'S2_rangesTested', 'ranges tested, before and after conversion to model units', ...
  'S3_fastingBenchmarks', 'Table S3: fasting benchmarks as reported and in mmol/day, with provenance notes', ...
  'syntheticCohort_uM', 'synthetic blood profiles (25 per sex, 26 metabolites, uM)', ...
  'S4_alanineThresholds', 'Table S4: alanine concentrations at which the hepatic alanine bounds cross', ...
  'S5_bmrByOrgan', 'Table S5: BMR change with heart rate and haematocrit by organ (kcal/day) and whole-body read-outs', ...
  'S6_heightCases', 'Table S6: height through cardiac output and blood volume', ...
  'S7_fastedSetupTests', 'Table S7: brain and heart ATP supply of the fasted models', ...
  'S8_fedFastedFluxes', 'Table S8: fed and fasted fluxes of the 11 read-outs (positive = release)', ...
  'S9_urineCreatinineFasted', 'Table S9: urine creatinine in the fasted models', ...
  'oatEffects', 'effect of every input on every read-out (FVA), with the shifts in minimum and maximum', ...
  'oatRanking', 'inputs ranked by their largest FVA effect', ...
  'oatL2Ranking', 'inputs ranked by reactions changed in the basal flux distribution, with BMR change', ...
  'feasibility', 'feasibility of every personalised model at each allowance', ...
  'crossedBounds', 'every crossed bound: reaction, metabolite, personalised and fixed value', ...
  'measurementPriority', 'feasibility after giving the measurement priority', ...
  'uptakeLimitRatio', 'median personalised / default organ uptake limit per metabolite', ...
  'activeState', 'activity capacity, energy expenditure, RQ and fuel exchanges per scenario', ...
  'boundMap_Harvey', 'every bound set by a measured blood concentration (value = k x concentration), Harvey', ...
  'boundMap_Harvetta', 'every bound set by a measured blood concentration, Harvetta');
writeBook(fullfile(outPath, 'SupplementaryData_Ch5_Inputs.xlsx'), X.in, desc);
writeBook(fullfile(outPath, 'SupplementaryData_Ch5_Results.xlsx'), X.res, desc);

tf = fopen(fullfile(P.latex, 'appendix_chapter5.tex'), 'w', 'n', 'UTF-8');
fprintf(tf, '%% Supplementary tables of Chapter 5, written by supplementaryMaterial.m\n');
fprintf(tf, '\\section*{Supplementary material}\n\\addcontentsline{toc}{section}{Supplementary material}\n');
fprintf(tf, '\\setcounter{table}{0}\n\\renewcommand{\\thetable}{S\\arabic{table}}\n\\renewcommand{\\theHtable}{S\\arabic{table}}\n\n');
fprintf(tf, ['The full inputs and results, including the synthetic profiles, every effect of every input on every read-out, ' ...
    'feasibility and crossed bounds of every personalised model and the bounds set by each measured concentration, are provided as Excel workbooks.\n\n']);
fprintf(tf, '%s\n', tex{:});
fclose(tf);
fprintf('Supplementary material written: %s and %s\n', fullfile(P.latex, 'appendix_chapter5.tex'), outPath);


%% helpers -------------------------------------------------------------------------------
function [st, rep] = reportedStat(stat, v1, v2, v3)
switch stat
    case 'meanSD',     st = 'mean $\pm$ SD';                  rep = sprintf('%g $\\pm$ %g', v1, v2);
    case 'meanCI',     st = 'mean (95\% CI)';                 rep = sprintf('%g (%g--%g)', v1, v2, v3);
    case 'pctl',       st = '2.5th--97.5th percentile';       rep = sprintf('%g--%g', v1, v2);
    case 'range',      st = '2.5th--97.5th percentile';       rep = sprintf('%g--%g', v1, v2);
    case 'meanP5P95',  st = 'mean (5th--95th percentile)';    rep = sprintf('%g (%g--%g)', v1, v2, v3);
    case 'lognP10P90', st = '10th, 50th, 90th percentile';    rep = sprintf('%g, %g, %g', v1, v2, v3);
    otherwise,         st = stat;                             rep = sprintf('%g, %g, %g', v1, v2, v3);
end
end

function u = unitTex(u)
map = {'umol/L', '$\mu$M'; 'mmol/L', 'mM'; 'mg/dL', 'mg\,dL$^{-1}$'; 'kg/m2', 'kg\,m$^{-2}$'; ...
       'mL/min/1.73m2', 'mL\,min$^{-1}$ per 1.73~m$^2$'; 'beats/min', 'beats\,min$^{-1}$'; 'L/L', ''; 'cm', 'cm'; ...
       'umol/kg/min', '$\mu$mol\,kg$^{-1}$\,min$^{-1}$'; 'mg/min/65kg', 'mg\,min$^{-1}$ per 65~kg'; 'mmol/day', 'mmol\,day$^{-1}$'};
i = find(strcmp(map(:, 1), u), 1);
if ~isempty(i), u = map{i, 2}; end
end

function v = S1fit(panel, refV, r, which)
q = strcmp(panel.vmh, refV.vmh{r});
if strcmp(refV.sex{r}, 'female'), p = [panel.p025Female(q) panel.p975Female(q)]; else, p = [panel.p025Male(q) panel.p975Male(q)]; end
v = p(which);
end

function s = fmt(x)
if ~isfinite(x), s = '--';
elseif abs(x) >= 100, s = sprintf('%.0f', x);
elseif abs(x) >= 10,  s = sprintf('%.1f', x);
elseif abs(x) >= 1,   s = sprintf('%.2f', x);
else,                 s = sprintf('%.2g', x);
end
s = strrep(s, '-', '$-$');
end

function s = signed(x)
% signed integer for a LaTeX table: +12, $-$12, 0
x = round(x);
if ~isfinite(x), s = '--'; elseif x == 0, s = '0'; elseif x > 0, s = sprintf('+%d', x); else, s = sprintf('$-$%d', -x); end
end

function s = ifNaN(x, txt)
if isnan(x), s = txt; else, s = sprintf('%.0f', x); end
end

function s = ifelseStr(c, a, b)
if c, s = a; else, s = b; end
end

function n = organName(n)
map = {'Adipocytes', 'Adipose tissue'; 'sIEC', 'Small intestine'; 'RBC', 'Red blood cells'; 'Agland', 'Adrenal gland'; ...
       'Scord', 'Spinal cord'; 'Gall', 'Gall bladder'; 'Thyroidgland', 'Thyroid gland'; 'Urinarybladder', 'Urinary bladder'};
i = find(strcmp(map(:, 1), n), 1);
if ~isempty(i), n = map{i, 2}; end
end

function s = latexTable(head, caption, label, hdr, align, rows)
% one float: \small table with booktabs rules, bold header and a bold caption heading
b = sprintf('\\begin{table}[H]\n\\centering\n\\small\n\\setlength{\\tabcolsep}{4pt}\n\\begin{tabular}{@{}%s@{}}\n\\toprule\n', align);
b = [b strjoin(strcat('\textbf{', hdr, '}'), ' & ') ' \\' newline '\midrule' newline];
for r = 1:size(rows, 1), b = [b strjoin(rows(r, :), ' & ') ' \\' newline]; end %#ok<AGROW>
b = [b '\bottomrule' newline '\end{tabular}' newline];
b = [b sprintf('\\caption[%s]{\\textbf{%s.} %s}\n\\label{%s}\n\\end{table}\n', head, head, caption, label)];
s = b;
end

function writeBook(file, S, desc)
if exist(file, 'file'), delete(file); end
names = fieldnames(S);
d = cellfun(@(n) desc.(n), names, 'UniformOutput', false);
writetable(table(string(names), string(d), 'VariableNames', {'sheet', 'contents'}), file, 'Sheet', 'README');
for i = 1:numel(names)
    T = S.(names{i});
    if ~istable(T), continue, end
    writetable(T, file, 'Sheet', names{i}(1:min(end, 31)));
end
end
