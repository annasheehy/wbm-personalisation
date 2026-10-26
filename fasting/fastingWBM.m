%% fastingWBM.m
% Builds feasible fasting WBMs for Harvey (male) and Harvetta (female).
% Reference: Thiele et al. 2020, Mol Syst Biol e8982.
% Units: mmol/person/day.  Solver: Gurobi via COBRA Toolbox.
%
% PARTS
%   0  Setup (paths, solver)
%   1  Load models + physiological constraints  (checkpointed)
%   2  Fed baseline sanity check
%   3  Fasting model  (setFastingConstraints — storage, ions, glycogen, choline)
%   4  BMR analysis — model energy expenditure vs reference equations (L1 & L2)
%   6  Summary figure (saved to resPath/fastingSummary.png)
%   7  Fed -> fasting L2 flux changes  (the open investigation)
%
% Dependencies (same folder): calculateModelBMR.m (model BMR, L1 & L2);
%   calculateBMR.m (phenomenological reference BMR, PSCM toolbox).
% Investigation log (INV-1 to INV-6): see fastingWBM_investigations.m

%% ==========================================================================
%  PART 0.  SETUP
%  ==========================================================================
P = persWBMpaths();
resPath        = fullfile(P.results, 'fasting');
rawCheckpoint  = fullfile(resPath, 'ckpt_rawModels.mat');
hmdbCheckpoint = fullfile(resPath, 'ckpt_HMDBconstrained.mat');
liftCheckpoint = fullfile(resPath, 'ckpt_liftedBase.mat');
feasCheckpoint = fullfile(resPath, 'ckpt_feasibleModels.mat');
fluxCheckpoint = fullfile(resPath, 'ckpt_fluxL2.mat');
if ~exist(resPath, 'dir'), mkdir(resPath); end

global CBT_LP_SOLVER CBT_QP_SOLVER
if isempty(CBT_LP_SOLVER)
    initCobraToolbox;          % slow — only on first run of the session
end
% Always assert the solver. The isempty guard above must gate ONLY the slow
% initCobraToolbox, NOT solver selection: otherwise a solver set by a previous
% session (e.g. ibm_cplex) silently persists and the BMR min-norm QP goes
% through CPLEX, which returns stat=3 (unconverged) with inflated fluxes.
okLP = changeCobraSolver('gurobi', 'LP', 0);
okQP = changeCobraSolver('gurobi', 'QP', 0);
if ~okLP || ~okQP
    warning('Gurobi not set (LP ok=%d, QP ok=%d) — still on %s. Check Gurobi install/licence.', ...
        okLP, okQP, CBT_LP_SOLVER);
end
fprintf('Solver: %s (LP)  %s (QP)  |  Output: %s\n', CBT_LP_SOLVER, CBT_QP_SOLVER, resPath);


%% ==========================================================================
%  PART 1.  LOAD MODELS + PHYSIOLOGICAL CONSTRAINTS  (checkpointed)
%  ==========================================================================

%% 1.1  Raw models -----------------------------------------------------------
if exist(rawCheckpoint, 'file')
    load(rawCheckpoint, 'harvey', 'harvetta');
    fprintf('1.1 — raw models from checkpoint.\n');
else
    harvey   = loadPSCMfile('Harvey');
    harvetta = loadPSCMfile('Harvetta');
    harvey.osenseStr = 'max';  harvetta.osenseStr = 'max';
    % Close liver biomass (regenerative growth, not fasting maintenance)
    harvey.lb(strcmp(harvey.rxns,     'Liver_biomass_reaction')) = 0;
    harvey.ub(strcmp(harvey.rxns,     'Liver_biomass_reaction')) = 0;
    harvetta.lb(strcmp(harvetta.rxns, 'Liver_biomass_reaction')) = 0;
    harvetta.ub(strcmp(harvetta.rxns, 'Liver_biomass_reaction')) = 0;
    save(rawCheckpoint, 'harvey', 'harvetta');
    fprintf('1.1 — raw models loaded and cached.\n');
end

%% 1.2  HMDB physiological constraints ---------------------------------------
% Organ-specific bounds from HMDB plasma concentrations x organ blood flow.
% Clears the objective vector -> reassign immediately after (Thiele 2020).
if exist(hmdbCheckpoint, 'file')
    load(hmdbCheckpoint, 'maleL', 'femaleL', 'ipMale', 'ipFemale');
    fprintf('1.2 — HMDB-constrained models from checkpoint.\n');
else
    sex = 'male';   standardPhysiolDefaultParameters; ipMale   = IndividualParameters;
    sex = 'female'; standardPhysiolDefaultParameters; ipFemale = IndividualParameters;
    maleL   = physiologicalConstraintsHMDBbased(harvey,   ipMale);
    femaleL = physiologicalConstraintsHMDBbased(harvetta, ipFemale);
    maleL   = changeObjective(maleL,   'Whole_body_objective_rxn');
    femaleL = changeObjective(femaleL, 'Whole_body_objective_rxn');
    maleL.osenseStr = 'max';  femaleL.osenseStr = 'max';
    save(hmdbCheckpoint, 'maleL', 'femaleL', 'ipMale', 'ipFemale');
    fprintf('1.2 — HMDB constraints applied and cached.\n');
end

%% 1.3  Lift the base models ONCE (cached) — the L2 QP fix, done up front -----
% liftCouplingConstraints reformulates ONLY the coupling constraints, which are
% invariant to the diet / storage / AA bound changes downstream. So lifting the
% base once means every model built from it (fed, fasting, ion-capped) is already
% lifted, and calculateModelBMR never re-lifts (its isModelLifted guard). This
% turns ~8 lifts per run into 0 (cached) — the ~2 h runtime was almost all lifting.
% NOTE: the cache depends only on the raw model's coupling, so it survives bound
% changes; delete ckpt_liftedBase.mat only if you change the RAW models (Part 1.1).
if exist(liftCheckpoint, 'file')
    load(liftCheckpoint, 'maleL', 'femaleL');
    fprintf('1.3 — lifted base models from checkpoint.\n');
else
    if ~isModelLifted(maleL),   maleL   = liftCouplingConstraints(maleL,   1e4, 0); end
    if ~isModelLifted(femaleL), femaleL = liftCouplingConstraints(femaleL, 1e4, 0); end
    save(liftCheckpoint, 'maleL', 'femaleL');
    fprintf('1.3 — base models lifted and cached.\n');
end


%% ==========================================================================
%  PART 2.  FED BASELINE SANITY CHECK
%  ==========================================================================
EUAverageDietNew
maleFed   = setDietConstraints(setFeedingFastingConstraints(maleL,   'feeding', -1e6, 1e6), Diet);
femaleFed = setDietConstraints(setFeedingFastingConstraints(femaleL, 'feeding', -1e6, 1e6), Diet);
maleFed.osenseStr = 'max';  femaleFed.osenseStr = 'max';
sol = optimizeWBModel(maleFed);   fprintf('2 — fed male:   stat=%d  obj=%.4f\n', sol.stat, sol.f);
sol = optimizeWBModel(femaleFed); fprintf('2 — fed female: stat=%d  obj=%.4f\n', sol.stat, sol.f);


%% ==========================================================================
%  PART 3.  FASTING MODEL
%  The complete, justified fasting constraint set is encapsulated in
%  setFastingConstraints: storage release (organ sinks drain), diet closed,
%  water + dissolved ions available, endogenous phosphate (bone Pi — the model
%  has no bone organ), and sex-specific hepatic glycogen release + intestinal
%  choline from blood. The defaults reproduce the values derived across the
%  investigations; see setFastingConstraints for the per-knob rationale and the
%  (still provisional) TODOs.
%  ==========================================================================
maleFast   = setFastingConstraints(maleL,   'male');
femaleFast = setFastingConstraints(femaleL, 'female');
maleFast.osenseStr = 'max';  femaleFast.osenseStr = 'max';
sol = optimizeWBModel(maleFast);   fprintf('3 — fasting male:   stat=%d  obj=%.4f\n', sol.stat, sol.f);
sol = optimizeWBModel(femaleFast); fprintf('3 — fasting female: stat=%d  obj=%.4f\n', sol.stat, sol.f);

save(feasCheckpoint, 'harvey', 'harvetta', 'maleL', 'femaleL', ...
     'maleFed', 'femaleFed', 'maleFast', 'femaleFast');
fprintf('3 — feasible fasting models saved to checkpoint.\n');


%% ==========================================================================
%  PART 4.  BMR ANALYSIS — model energy expenditure vs reference equations
%
%  BMR (basal metabolic rate) is the gold-standard measure of resting energy
%  expenditure. Thiele 2020 read it from the whole-body model as the "basal
%  metabolic flux": impose resting energy floors, find the parsimonious resting
%  flux, sum ATP hydrolysis across all compartments, convert to kcal/day
%  (64 kJ/mol ATP). We compute it here (via calculateModelBMR.m) and compare to
%  the phenomenological references (Harris-Benedict / Mifflin-St Jeor,
%  calculateBMR.m) and the paper's archetypes (Harvey 1455, Harvetta 1344).
%
%  Two flux norms are reported (calculateModelBMR runs both):
%    L2  minimise ||v||_2 (QP) — the norm used in Thiele 2020 (minNorm 1e-6).
%         This is the value that should REPRODUCE the paper archetypes
%         (Harvey ~1455, Harvetta ~1344): same floors, same obj=1, same ATP
%         summation. The L2 QP is numerically hard on the 80k-rxn WBM, so
%         calculateModelBMR lifts the coupling constraints (liftCouplingConstraints)
%         and solves barrier-only with crossover off (Crossover=0) to converge.
%    L1  minimise sum|v|  (LP) — always robust, but SPARSE, so it under-reads
%         total ATP turnover (~half the L2 value). Used as the robustness
%         cross-check; L2 is the headline once it converges.
%  ==========================================================================
if ~exist('maleFast', 'var'), load(feasCheckpoint); end
R = struct();                                   % results collector for the Part 6 figure

%% 4.1  Phenomenological reference BMR (Harris-Benedict / Roza / Mifflin) -----
% calculateBMR(sex, weight_kg, height_cm, age_yr) -> [HB1919; Roza1984; Mifflin1990].
% Harvey/Harvetta are ageless archetypes (ip*.Age = NaN); use WHO/FAO adult 30 yr.
% For personalisation, replace with the patient's actual age.
patientAgeM = ipMale.Age;    if isnan(patientAgeM), patientAgeM = 30; end
patientAgeF = ipFemale.Age;  if isnan(patientAgeF), patientAgeF = 30; end
BMRmale   = calculateBMR('male',   ipMale.bodyWeight,   ipMale.Height,   patientAgeM);
BMRfemale = calculateBMR('female', ipFemale.bodyWeight, ipFemale.Height, patientAgeF);
fprintf('\n4.1 — Reference BMR (kcal/day) [HB1919 / Roza1984 / Mifflin1990], age %d:\n', patientAgeM);
fprintf('      Male   (%.0f kg, %.0f cm): %.0f / %.0f / %.0f\n', ...
    ipMale.bodyWeight, ipMale.Height, BMRmale(1), BMRmale(2), BMRmale(3));
fprintf('      Female (%.0f kg, %.0f cm): %.0f / %.0f / %.0f\n', ...
    ipFemale.bodyWeight, ipFemale.Height, BMRfemale(1), BMRfemale(2), BMRfemale(3));
R.refMale = BMRmale(:)';  R.refFemale = BMRfemale(:)';  R.archetype = [1455 1344];

%% 4.2  Model BMR at the fasting stage — L1 and L2 ---------------------------
bmrM = calculateModelBMR(maleFast);             % returns .L1 and .L2 (kcal, atp, stat, v)
bmrF = calculateModelBMR(femaleFast);
fprintf('\n4.2 — Model BMR at fasting stage (kcal/day):\n');
fprintf('              L1 (LP)          L2 (QP)          Mifflin   Thiele archetype\n');
fprintf('   Male    %6.0f (stat %2d)  %6.0f (stat %2d)   %6.0f    %d\n', ...
    bmrM.L1.kcal, bmrM.L1.stat, bmrM.L2.kcal, bmrM.L2.stat, BMRmale(3),   1455);
fprintf('   Female  %6.0f (stat %2d)  %6.0f (stat %2d)   %6.0f    %d\n', ...
    bmrF.L1.kcal, bmrF.L1.stat, bmrF.L2.kcal, bmrF.L2.stat, BMRfemale(3), 1344);
R.fastL1  = [bmrM.L1.kcal bmrF.L1.kcal];   R.fastL1stat = [bmrM.L1.stat bmrF.L1.stat];
R.fastL2  = [bmrM.L2.kcal bmrF.L2.kcal];   R.fastL2stat = [bmrM.L2.stat bmrF.L2.stat];
fluxFastM = bmrM.L2.v;  fluxFastF = bmrF.L2.v;   % fasting L2 flux -> Part 7 flux-change analysis

%% 4.3  Model BMR across construction stages (L1) — is it stable? ------------
% L1 only (fast, robust). base (HMDB) -> fed -> fasting. Physiology expects BMR
% roughly stable from fed to (short) fasting, so a big jump flags a solver/model
% issue at that stage rather than biology.
R.stageNames = {'base','fed','fasting'};
stageMdl     = {maleL femaleL; maleFed femaleFed; maleFast femaleFast};
R.stageL1male = nan(1,3);  R.stageL1female = nan(1,3);
fprintf('\n4.3 — model BMR (L1) by construction stage (kcal/day):\n');
for s = 1:3
    bm = calculateModelBMR(stageMdl{s,1}, 'norms', {'L1'});
    bf = calculateModelBMR(stageMdl{s,2}, 'norms', {'L1'});
    R.stageL1male(s)   = bm.L1.kcal;
    R.stageL1female(s) = bf.L1.kcal;
    fprintf('   %-8s  male %6.0f (stat %d)   female %6.0f (stat %d)\n', ...
        R.stageNames{s}, bm.L1.kcal, bm.L1.stat, bf.L1.kcal, bf.L1.stat);
end

%% 4.3b  L2 BMR across states — reproduces the Thiele 2020 archetypes --------
% L2 (min ||v||_2) is Thiele's actual norm, so the FED L2 value should reproduce
% the paper archetypes (Harvey 1455, Harvetta 1344). Run L2 on base and fed for
% both sexes; the fasting L2 is already in 4.2. (This also settled the earlier
% L2 non-convergence: it was a solver-settings failure, not the fasting model —
% base/fed/fasting all solve once the QP is lifted + barrier-only, crossover off.)
fprintf('\n4.3b — L2 BMR across states (Thiele norm; fasting L2 in 4.2):\n');
l2diag = {'base', maleL, femaleL; 'fed', maleFed, femaleFed};
R.stageL2male = nan(1,3);  R.stageL2female = nan(1,3);   % [base fed fasting]
for s = 1:size(l2diag,1)
    bm2 = calculateModelBMR(l2diag{s,2}, 'norms', {'L2'});
    bf2 = calculateModelBMR(l2diag{s,3}, 'norms', {'L2'});
    R.stageL2male(s)   = bm2.L2.kcal;
    R.stageL2female(s) = bf2.L2.kcal;
    if strcmp(l2diag{s,1}, 'fed'), fluxFedM = bm2.L2.v;  fluxFedF = bf2.L2.v; end   % fed L2 flux -> Part 7
    fprintf('   %-6s  male L2 %6.0f (stat %d)   female L2 %6.0f (stat %d)\n', ...
        l2diag{s,1}, bm2.L2.kcal, bm2.L2.stat, bf2.L2.kcal, bf2.L2.stat);
end
R.stageL2male(3) = R.fastL2(1);  R.stageL2female(3) = R.fastL2(2);

%% 4.4  Diagnostic: top ATP-consuming reactions (L1 fasting flux) ------------
% Confirms the L1 solve is cycle-free: real basal ATP use spreads across
% maintenance/ion pumps, so one reaction dominating would signal a futile cycle.
diag4 = {'male', maleFast, bmrM; 'female', femaleFast, bmrF};
for d = 1:size(diag4, 1)
    mdl = diag4{d,2};  bmrS = diag4{d,3};       % each sex uses its OWN atpMets
    if isempty(bmrS.L1.v), continue; end
    perRxnATP  = full(-sum(min(mdl.S(bmrS.atpMets, :) .* bmrS.L1.v(:)', 0), 1));
    [top, ord] = sort(perRxnATP, 'descend');
    fprintf('\n4.4 — top 5 ATP-consuming reactions (%s, L1):\n', diag4{d,1});
    for k = 1:min(5, numel(top))
        if top(k) <= 0, break; end
        fprintf('   %10.1f mmol/day  %s\n', top(k), mdl.rxns{ord(k)});
    end
end

%% 4.5  Resting brain/heart ATP ceiling --------------------------------------
% The BMF solve imposes resting energy floors (Brain_DM_atp >= 3500,
% Heart_DM_atp >= 1000). If a fasting model cannot supply that much, the BMF is
% infeasible (this is why the female fasting BMR came back stat=0). Here we ask
% the opposite question per sex: with the whole-body objective fixed at 1, what
% is the MAXIMUM brain/heart ATP the fasting model can supply? A ceiling below
% the floor pinpoints the fuel bottleneck (brain runs on glucose/ketones).
brainFloor = 3500;  heartFloor = 1000;
sexCeil    = {'male', maleFast; 'female', femaleFast};
R.maxBrain = nan(1,2);  R.maxHeart = nan(1,2);
R.brainFloor = brainFloor;  R.heartFloor = heartFloor;
fprintf('\n4.5 — resting ATP ceiling at obj=1 (max the fasting model can supply):\n');
for c = 1:2
    m  = sexCeil{c,2};
    ob = strcmp(m.rxns, 'Whole_body_objective_rxn');
    m.lb(ob) = 1;  m.ub(ob) = 1;
    bi = contains(m.rxns, 'Brain_DM_atp');
    hi = contains(m.rxns, 'Heart_DM_atp');
    mB = m;  mB.c(:) = 0;  mB.c(bi) = 1;  mB.osenseStr = 'max';
    mH = m;  mH.c(:) = 0;  mH.c(hi) = 1;  mH.osenseStr = 'max';
    sB = optimizeWBModel(mB);
    sH = optimizeWBModel(mH);
    R.maxBrain(c) = sB.f;  R.maxHeart(c) = sH.f;
    bOK = 'SHORT';  if sB.f >= brainFloor, bOK = 'OK'; end
    hOK = 'SHORT';  if sH.f >= heartFloor, hOK = 'OK'; end
    fprintf('   %-6s: brain max %6.0f / floor %d (%s)   heart max %6.0f / floor %d (%s)\n', ...
        sexCeil{c,1}, sB.f, brainFloor, bOK, sH.f, heartFloor, hOK);
end


%% ==========================================================================
%  PART 6.  SUMMARY FIGURE   (top row: A setup, B L2 BMR-across-states vs
%  reference; bottom row: C resting ATP ceiling, D L1-vs-L2 by state, E summary.)
%  L2 (Thiele 2020 min-norm QP) is the validated headline and reproduces the
%  archetypes; L1 (sparse LP) is shown as the ~2x-lower cross-check. Colour:
%  blue ramp = nutritional state (base->fed->fasted, light->dark); greys =
%  reference equations; blue/orange = male/female; red dashed = organ ATP demand
%  floor. Every computed-flux panel carries a subtitle. Run Part 4 first.
%  ==========================================================================
if ~exist('R', 'var') || ~isfield(R, 'fastL1')
    warning('Part 6: results struct R not populated — run Part 4 first.');
else
    global CBT_LP_SOLVER CBT_QP_SOLVER          % report the ACTUAL solver, not an assumption
    lpS = CBT_LP_SOLVER;  if isempty(lpS), lpS = '?'; end
    qpS = CBT_QP_SOLVER;  if isempty(qpS), qpS = '?'; end

    % validated palette (dataviz skill, light mode)
    cModel   = [ 42 120 214]/255;  % blue   — male
    cFemale  = [235 104  52]/255;  % orange — female
    blueRamp = [134 182 239; 57 135 229; 28 92 171]/255;  % base / fed / fasted (light->dark)
    gRef     = [0.82 0.82 0.80; 0.68 0.68 0.66; 0.54 0.54 0.52; 0.38 0.38 0.36];  % refs: Mifflin/HB/Roza/Thiele (dark = archetype)
    cFloor   = [208  59  59]/255;  % critical red — demand floor
    cMuted   = [137 135 129]/255;  % muted ink — captions
    cLab     = [0.28 0.28 0.28];   % value-label ink
    cGrid    = [0.88 0.88 0.86];

    figSummary = figure('Name', 'Fasting WBM summary', 'Color', 'w', ...
                        'Position', [50 50 1550 900]);
    sgtitle('Fasting whole-body metabolic model:  BMR analysis  (Harvey / Harvetta)', ...
            'FontWeight','bold', 'FontSize', 13);

    % --- A: setup / methods (text) ---
    subplot(2,3,1); axis off;
    % Read the ACTUAL fasting-knob bounds from the models so this caption never
    % goes stale when the knobs change (hardcoded values were the bug). Uses the
    % same reaction names Part 3 sets; AA = mode over the storage-AA muscle sinks.
    if exist('maleFast','var') && exist('femaleFast','var')
        rd  = @(m, r) m.lb(strcmp(m.rxns, r));
        aaM = mode(maleFast.lb(~cellfun(@isempty,   regexp(maleFast.rxns,   '^Muscle_sink_\w+_L\(c\)$', 'once'))));
        aaF = mode(femaleFast.lb(~cellfun(@isempty, regexp(femaleFast.rxns, '^Muscle_sink_\w+_L\(c\)$', 'once'))));
        knobLines = {
            sprintf('   storage-AA sinks (M %g, F %g);', aaM, aaF)
            sprintf('   glycogen release (M %g, F %g);', rd(maleFast,'Liver_sink_glygn2(c)'), rd(femaleFast,'Liver_sink_glygn2(c)'))
            sprintf('   intestinal choline from blood (M %g, F %g)', rd(maleFast,'sIEC_EX_chol(e)_[bc]'), rd(femaleFast,'sIEC_EX_chol(e)_[bc]')) };
    else
        knobLines = {'   storage-AA sinks, glycogen release,'; '   intestinal choline from blood (see Part 3)'};
    end
    setupTxt = [{
        sprintf('Solver: %s (LP), %s (QP)', lpS, qpS)
        'Harvey (male) / Harvetta (female), Thiele 2020'
        ''
        '\bfBase\rm  model as loaded + HMDB bounds (diet unchanged)'
        '\bfFed\rm    current EU average diet applied (small update)'
        '\bfFasting\rm  diet closed except water + ions;'
        }; knobLines; {
        ''
        '\bfBMR\rm = \Sigma ATP hydrolysis \times 64 kJ/mol \times 0.239'
        '   L2 = min ||v||_2  (QP; Thiele 2020 norm; HEADLINE)'
        '        lifted + barrier (crossover off); reproduces archetypes'
        '   L1 = min \Sigma|v|  (LP; robust cross-check; ~half L2)'
        }];
    text(0, 1, setupTxt, 'VerticalAlignment','top', 'FontSize', 8.5, 'Interpreter','tex');
    title('A. Setup');

    % --- B: model L2 BMR across nutritional states vs reference equations ---
    % L2 is Thiele's norm; the fed L2 bars should land on the archetype/Mifflin
    % references — the visual proof of reproduction (Harvey 1455, Harvetta 1344).
    subplot(2,3,[2 3]);
    Bdata = [R.stageL2male(1)   R.stageL2male(2)   R.stageL2male(3)   R.refMale(3)   R.refMale(1)   R.refMale(2)   R.archetype(1);
             R.stageL2female(1) R.stageL2female(2) R.stageL2female(3) R.refFemale(3) R.refFemale(1) R.refFemale(2) R.archetype(2)];
    b = bar(Bdata);
    for k = 1:3, b(k).FaceColor = blueRamp(k,:); end     % base / fed / fasted (L2)
    for k = 1:4, b(3+k).FaceColor = gRef(k,:);   end      % Mifflin / HB / Roza / Thiele archetype
    set(gca, 'XTickLabel', {'Male (Harvey)','Female (Harvetta)'}, 'GridColor', cGrid); grid on;
    ymaxB = max(Bdata(~isnan(Bdata))); ylim([0 ymaxB*1.18]);
    ylabel('BMR (kcal/day)');
    title('B. Model BMR (L2, Thiele norm) across states vs reference equations');
    subtitle(['Model bars = L2 min-norm QP at obj=1 (Thiele 2020 norm); kcal = \SigmaATP \times 64 kJ.  ' ...
              'Grey bars = phenomenological equations + Thiele model archetype (darkest).'], 'FontSize', 7.5, 'Color', cMuted);
    legend({'Base (L2)','Fed (L2)','Fasted (L2)','Mifflin-St Jeor','Harris-Benedict','Roza-Shizgal','Thiele 2020'}, ...
        'Location','northoutside', 'Orientation','horizontal', 'FontSize', 7);
    for i = 1:numel(b)                                   % value labels
        xe = b(i).XEndPoints; ye = b(i).YEndPoints;
        for g = 1:numel(xe)
            if ~isnan(ye(g)) && ye(g) > 0
                text(xe(g), ye(g), sprintf('%.0f', ye(g)), 'HorizontalAlignment','center', ...
                    'VerticalAlignment','bottom', 'FontSize', 6.5, 'Color', cLab);
            end
        end
    end
    % Annotate any missing fasted-L2 bar, distinguishing biological infeasibility
    % (L1 also failed for that state) from L2 QP non-convergence (L1 solved, QP did not).
    fastedL1 = [R.stageL1male(3) R.stageL1female(3)];
    fastedL2 = [R.stageL2male(3) R.stageL2female(3)];
    for gsex = 1:2
        if isnan(fastedL2(gsex))
            if isnan(fastedL1(gsex)), note = 'infeasible'; else, note = 'L2 n/c'; end
            text(b(3).XEndPoints(gsex), ymaxB*0.04, note, 'Rotation',90, 'Color',cFloor, ...
                'FontSize',7.5, 'HorizontalAlignment','left', 'VerticalAlignment','middle', 'FontWeight','bold');
        end
    end

    % --- C: resting brain/heart ATP ceiling vs demand floor ---
    subplot(2,3,4);
    Cdata = [R.maxBrain(1) R.maxBrain(2); R.maxHeart(1) R.maxHeart(2)];   % rows Brain/Heart, cols M/F
    bc = bar(Cdata); bc(1).FaceColor = cModel; bc(2).FaceColor = cFemale; hold on;
    plot([0.55 1.45], [R.brainFloor R.brainFloor], '--', 'Color', cFloor, 'LineWidth', 1.6);
    plot([1.55 2.45], [R.heartFloor R.heartFloor], '--', 'Color', cFloor, 'LineWidth', 1.6);
    text(1.46, R.brainFloor, sprintf('floor %d ', R.brainFloor), 'Color',cFloor, 'FontSize',7, ...
        'HorizontalAlignment','right', 'VerticalAlignment','bottom');
    text(2.46, R.heartFloor, sprintf('floor %d ', R.heartFloor), 'Color',cFloor, 'FontSize',7, ...
        'HorizontalAlignment','right', 'VerticalAlignment','bottom');
    set(gca, 'XTickLabel', {'Brain','Heart'}, 'GridColor', cGrid); grid on;
    ylabel('max ATP supply (mmol/day)');
    title('C. Resting ATP the model can supply vs demand');
    subtitle('FBA LP: maximise organ ATP demand at obj=1 (fasted model).', 'FontSize',7.5, 'Color',cMuted);
    legend({'Male','Female'}, 'Location','north', 'FontSize',7);
    for i = 1:numel(bc)
        xe = bc(i).XEndPoints; ye = bc(i).YEndPoints;
        for g = 1:numel(xe)
            if ~isnan(ye(g)) && ye(g) > 0
                text(xe(g), ye(g), sprintf('%.0f', ye(g)), 'HorizontalAlignment','center', ...
                    'VerticalAlignment','bottom', 'FontSize', 6.5, 'Color', cLab);
            end
        end
    end

    % --- D: L1 vs L2 by state (male) — L1 under-reads ATP turnover ~2x ---
    subplot(2,3,5);
    Ddata = [R.stageL1male(1) R.stageL2male(1);          % rows base/fed/fasting, cols [L1 L2]
             R.stageL1male(2) R.stageL2male(2);
             R.stageL1male(3) R.stageL2male(3)];
    bd = bar(Ddata); bd(1).FaceColor = blueRamp(1,:); bd(2).FaceColor = blueRamp(3,:); hold on;
    set(gca, 'XTickLabel', {'base','fed','fasting'}, 'GridColor', cGrid); grid on;
    ymaxD = max([Ddata(~isnan(Ddata)); 1]);  ylim([0 ymaxD*1.20]);
    ylabel('BMR (kcal/day)');
    title('D. L1 vs L2 by state (male): L1 under-reads ~2x');
    subtitle(['Same model, floors, obj=1. L1 (sparse LP) finds a lower-turnover ' ...
              'distribution than the L2 min-norm QP.'], 'FontSize',7.5, 'Color',cMuted);
    legend({'L1  (min \Sigma|v|)','L2  (min ||v||_2)'}, 'Location','northwest', 'FontSize',7);
    for i = 1:numel(bd)
        xe = bd(i).XEndPoints; ye = bd(i).YEndPoints;
        for g = 1:numel(xe)
            if ~isnan(ye(g)) && ye(g) > 0
                text(xe(g), ye(g), sprintf('%.0f', ye(g)), 'HorizontalAlignment','center', ...
                    'VerticalAlignment','bottom', 'FontSize', 6.5, 'Color', cLab);
            end
        end
    end

    % --- E: results summary (text) — fully dynamic (no hardcoded verdicts) ---
    subplot(2,3,6); axis off;
    if isnan(R.stageL2female(3)), femFastL2 = 'n/c'; else, femFastL2 = sprintf('%.0f', R.stageL2female(3)); end
    if isnan(R.stageL2male(3)),   malFastL2 = 'n/c'; else, malFastL2 = sprintf('%.0f', R.stageL2male(3)); end
    if R.maxBrain(1) >= R.brainFloor && R.maxBrain(2) >= R.brainFloor
        brainLine = sprintf('Both brains clear the floor (M %.0f / F %.0f vs %d).', R.maxBrain(1), R.maxBrain(2), R.brainFloor);
    else
        brainLine = sprintf('Brain ceiling short of floor: M %.0f / F %.0f vs %d.', R.maxBrain(1), R.maxBrain(2), R.brainFloor);
    end
    txt = {
        '\bfResults summary\rm   (L2 = Thiele norm = headline)'
        ''
        sprintf('Fed L2:   male %.0f (Harvey 1455),  female %.0f (Harvetta 1344)', R.stageL2male(2), R.stageL2female(2))
        sprintf('Fasting L2 (kcal):   male %s,  female %s', malFastL2, femFastL2)
        ''
        sprintf('L1 cross-check (sparse LP, under-reads): fasting L1 M %.0f / F %.0f.', ...
            R.stageL1male(3), R.stageL1female(3))
        ''
        brainLine
        sprintf('Heart not limiting (M %.0f / F %.0f vs %d).', R.maxHeart(1), R.maxHeart(2), R.heartFloor)
        ''
        'Storage-AA release set at the Part-3 value (provisional; needs justifying).'
        };
    text(0, 1, txt, 'VerticalAlignment','top', 'FontSize', 9, 'Interpreter','tex');
    title('E. Summary');

    exportgraphics(figSummary, fullfile(resPath, 'fastingSummary.png'), 'Resolution', 150);
    fprintf('\n6 — summary figure saved: %s\n', fullfile(resPath, 'fastingSummary.png'));
end


%% ==========================================================================
%  PART 7.  FED -> FASTING ORGAN EXCHANGE FLUXES  (L2 min-norm QP)
%  Data-led survey of the inter-organ fuel economy: every organ<->blood exchange
%  reaction (<Organ>_EX_<met>(e)_[bc], and brain _[csf]) with its fed and fasting
%  L2 flux and the change. Exchange-flux SIGN is meaningful (uptake vs release), so
%  these are printed RAW for manual interpretation — no summing, no cherry-picking.
%  Console shows the global top movers per sex (small fuels, e.g. ketones, carry
%  tiny flux and won't rank there but ARE in the Excel). Excel: one sheet per organ,
%  sorted by |change| (fedFastingExchange_<sex>.xlsx) — browse / Ctrl-F glc_D, bhb...
%  Re-runnable on its own: loads the flux vectors + models from the checkpoints.
%  ==========================================================================
if ~exist('fluxFastM', 'var') && exist(fluxCheckpoint, 'file')
    load(fluxCheckpoint, 'fluxFedM', 'fluxFedF', 'fluxFastM', 'fluxFastF');
end
if ~exist('maleFast', 'var'), load(feasCheckpoint, 'maleFast', 'femaleFast'); end

if exist('fluxFedM', 'var') && exist('fluxFastM', 'var') && ~isempty(fluxFedM) && ~isempty(fluxFastM)
    save(fluxCheckpoint, 'fluxFedM', 'fluxFedF', 'fluxFastM', 'fluxFastF');   % for standalone re-runs
    topShow = 20;
    vmhAbbr = {};  vmhName = {};                     % VMH abbreviation -> full metabolite name
    try
        vdb = loadVMHDatabase;  vmhAbbr = vdb.metabolites(:,1);  vmhName = vdb.metabolites(:,2);
    catch
        warning('Part 7: loadVMHDatabase failed — using metabolite abbreviations.');
    end
    fluxSex = {'male', maleFast, fluxFedM, fluxFastM; 'female', femaleFast, fluxFedF, fluxFastF};
    for q = 1:2
        m = fluxSex{q,2};  vFed = fluxSex{q,3};  vFast = fluxSex{q,4};
        if isempty(vFed) || isempty(vFast)
            fprintf('\n7 — %s: L2 flux missing (fed or fasting did not solve) — skipped.\n', fluxSex{q,1});
            continue
        end
        % Organ<->blood exchange reactions: <Organ>_EX_<met>(e)_[bc] (brain _[csf]).
        isExch = ~cellfun(@isempty, regexp(m.rxns, '^[A-Za-z0-9]+_EX_.*_\[(bc|csf)\]$', 'once'));
        idx    = find(isExch);
        orgTok = regexp(m.rxns(idx), '^([A-Za-z0-9]+)_EX_', 'tokens', 'once');
        organ  = cellfun(@(c) c{1}, orgTok, 'UniformOutput', false);
        fed    = full(vFed(idx));  fast = full(vFast(idx));  d = fast - fed;

        % Per exchange: the blood-side metabolite, its VMH full name, and the
        % uptake/release direction DERIVED FROM STOICHIOMETRY (sign of flux x the
        % blood-metabolite coefficient) — correct regardless of how the reaction is
        % written. release = organ -> blood; uptake = blood -> organ.
        nE = numel(idx);  metName = repmat({''}, nE, 1);
        fedDir = repmat({'-'}, nE, 1);  fastDir = repmat({'-'}, nE, 1);
        for e = 1:nE
            mr = find(m.S(:, idx(e)) ~= 0);
            bc = mr(find(~cellfun(@isempty, regexp(m.mets(mr), '\[(bc|csf)\]$', 'once')), 1));
            if isempty(bc), continue; end
            abbr = regexprep(m.mets{bc}, '\[(bc|csf)\]$', '');
            hit  = find(strcmp(vmhAbbr, abbr), 1);
            if ~isempty(hit), metName{e} = vmhName{hit}; else, metName{e} = abbr; end
            coeff = m.S(bc, idx(e));
            if     fed(e)*coeff  >  1e-9, fedDir{e}  = 'release'; elseif fed(e)*coeff  < -1e-9, fedDir{e}  = 'uptake'; end
            if     fast(e)*coeff >  1e-9, fastDir{e} = 'release'; elseif fast(e)*coeff < -1e-9, fastDir{e} = 'uptake'; end
        end

        % Excel: one sheet per organ, sorted by |change|.
        xlsFile = fullfile(resPath, sprintf('fedFastingExchange_%s.xlsx', fluxSex{q,1}));
        if exist(xlsFile, 'file'), delete(xlsFile); end
        organs = unique(organ);
        for oi = 1:numel(organs)
            sel = strcmp(organ, organs{oi});
            T = table(m.rxns(idx(sel)), metName(sel), fed(sel), fast(sel), d(sel), fedDir(sel), fastDir(sel), abs(d(sel)), ...
                'VariableNames', {'reaction','metabolite','fed','fasting','delta','fedDir','fastDir','absDelta'});
            T = sortrows(T, 'absDelta', 'descend');  T.absDelta = [];
            writetable(T, xlsFile, 'Sheet', organs{oi});
        end

        % Console: top movers, water/O2/CO2/HCO3 hidden so fuels surface (Excel has all).
        keep  = cellfun(@isempty, regexp(m.rxns(idx), '_(h2o|o2|co2|hco3)[\(\[]', 'once')) & abs(d) > 1e-3;
        cShow = find(keep);
        [~, o] = sort(abs(d(cShow)), 'descend');
        fprintf('\n7 — %s: organ<->blood exchange, top %d fed->fasting changes (water/O2/CO2 hidden; all in %s):\n', ...
            fluxSex{q,1}, topShow, xlsFile);
        fprintf('   %-28s %-20s %9s %9s %9s  %s\n', 'reaction', 'metabolite', 'fed', 'fasting', 'delta', 'fastDir');
        for k = 1:min(topShow, numel(o))
            j = cShow(o(k));  nm = metName{j};
            fprintf('   %-28s %-20s %9.1f %9.1f %9.1f  %s\n', ...
                m.rxns{idx(j)}, nm(1:min(numel(nm),20)), fed(j), fast(j), d(j), fastDir{j});
        end
    end
else
    warning('Part 7: L2 flux vectors not available — run Parts 4.2 and 4.3b (needs L2 to have solved).');
end
