%% fastingWBM_investigations.m
% Archive of feasibility investigations carried out while building the fasting WBM.
% Each section is standalone — it loads models from the checkpoint if needed.
% Results are documented in comments so sections do not need to be re-run.
%
% SUMMARY OF FINDINGS
%   INV-1  findViableDiet on infeasible fasting model (-10 sinks)
%          -> flags 7thf, ksi/ksii, choline (female). All are artefacts of -10 sink cap.
%   INV-2  7thf demand is really a glutamate shortage from the -10 sink cap
%   INV-3  ksi/ksii demand: missing cartilage organ in Harvey (raw model artefact)
%   INV-4  Female choline shortfall: sIEC uptake cap -14.93 too tight
%   INV-5  noTrTr objective: only 4/2 organs have it; cannot rewire whole-body objective
%   INV-6  Bare minimum: -1e6 sinks + ions + sIEC choline -> feasible (both sexes)
%          Threshold: female ~-20, male ~-30 (multiples of 10; exact values TBD)
%
% Current feasible model approach: see fastingWBM.m Part 3.


%% ==========================================================================
%  SETUP  (run Part 0 of fastingWBM.m first, or set these paths manually)
%  ==========================================================================
if ~exist('feasCheckpoint', 'var')
    P = persWBMpaths();
    resPath        = fullfile(P.results, 'fasting');
    feasCheckpoint = fullfile(resPath, 'ckpt_feasibleModels.mat');
end


%% ==========================================================================
%  INV-1.  findViableDiet on infeasible fasting model (-10 sink cap)
%
%  Question: what is the minimum dietary addition that restores feasibility
%  when sinks are at the -10 default?
%
%  RESULT (documented — slow to re-run):
%    Male:   Diet_EX_7thf (-39.6)   Diet_EX_ksi (-0.07)   Diet_EX_ksii_core4 (-13.1)
%    Female: Diet_EX_7thf (-24.8)   Diet_EX_chol (-1.27)  Diet_EX_ksi (-0.07)
%            Diet_EX_ksii_core4 (-11.2)
%  Three issues flagged: 7thf -> INV-2, ksi/ksii -> INV-3, chol (female) -> INV-4.
%  ==========================================================================
fprintf('\nINV-1 — documented result: findViableDiet flags 7thf (both), ksi/ksii (both), chol (female).\n');
% Uncomment to re-run (slow, ~5 min each):
% if ~exist('maleFastInfeasible','var'), load(feasCheckpoint); end
% [~,~,~,dvdM] = findViableDiet(maleFastInfeasible, 'Diet_EX_'); disp(dvdM)
% [~,~,~,dvdF] = findViableDiet(femaleFastInfeasible, 'Diet_EX_'); disp(dvdF)


%% ==========================================================================
%  INV-2.  7thf demand is a glutamate shortage, not a folate issue
%
%  7thf = heptaglutamyl-tetrahydrofolate: a retained folate carrying 7 glutamate
%  residues. Providing THF (folate) does NOT remove the demand; providing glutamate
%  (280 mmol = 40 molecules x 7 residues) removes it completely in both sexes.
%  The model pulls 7thf purely for its glutamate — the folate is incidental cargo.
%
%  ROOT CAUSE: setFeedingFastingConstraints caps all 20 muscle AA sinks at -10,
%  which is far too tight for glutamate (the most abundant free AA in muscle).
%
%  RESULT: binary search found minimum Muscle_sink_glu_L(c) for male feasibility
%  is ~-221 mmol/day. Felig 1970 fasting muscle glutamate release: ~60-150 mmol/day.
%  => model demands more than physiology — likely a missing interorgan route.
%
%  OPEN QUESTION: female-specific minimum not yet determined (-280 is a placeholder).
%  ==========================================================================
if ~exist('maleFastInfeasible','var'), load(feasCheckpoint); end

if ~exist('gluMin', 'var')
    % Binary search: male, with GAG proxy and choline bound already fixed
    gluBase = changeRxnBounds(maleFastInfeasible, 'sIEC_EX_chol(e)_[bc]', -20, 'l');
    for g = {'ksi','ksii_core2','ksii_core4'}
        gluBase = addReaction(gluBase, ['Cartilage_EX_' g{1} '[bc]'], ...
            'reactionFormula', ['-> ' g{1} '[bc]'], 'lowerBound', 0, 'upperBound', 1e6);
    end
    lo = -280; hi = -10;
    for it = 1:12
        mid = (lo + hi)/2;
        s = optimizeWBModel(changeRxnBounds(gluBase, 'Muscle_sink_glu_L(c)', mid, 'l'));
        if s.stat == 1, lo = mid; else, hi = mid; end
    end
    gluMin = lo;
end
fprintf('\nINV-2 — min Muscle_sink_glu_L(c) for male feasibility: %.1f mmol/day\n', gluMin);
fprintf('        Felig 1970 physiological range: 60-150 mmol/day\n');
fprintf('        Model over-demands glutamate — missing interorgan pathway suspected.\n');


%% ==========================================================================
%  INV-3.  ksi / ksii keratan sulfate: missing cartilage organ artefact
%
%  ksi / ksii_core2 / ksii_core4 are cartilage glycosaminoglycans. Harvey and
%  Harvetta have NO bone/cartilage organ (Thiele 2020 known limitation). Colon
%  and Liver carry minimum uptake bounds for ksii_core4 from the RAW model
%  construction (not HMDB). With no cartilage organ there is no source in fasting.
%
%  RESULT — origin trace:
%    Colon_EX_ksii_core4:  lb = -4.4554 at raw, post-HMDB, and post-fasting (unchanged)
%    Liver_EX_ksii_core4:  lb = -9.2621 at raw, post-HMDB, and post-fasting (unchanged)
%    Kidney_EX_ksii_core4: lb = -3.7368 raw -> -2.5920 post-HMDB (relaxed by HMDB)
%  => Colon/Liver bounds come from the RAW model, not HMDB.
%
%  FVA of proxy reactions: ksi 0-0.004, ksii_core2 9.13-9.46, ksii_core4 14.57-14.84
%  Published cartilage GAG turnover: ~0.01-0.05 mmol/day => model is 300-1500x too high.
%
%  OPEN QUESTION: Harvetta is feasible WITHOUT the proxy — why does only Harvey need it?
%  This was the fix applied in the OLD script (now removed). Under the new approach
%  (-30/-20 sink bounds) the proxy reactions are NOT needed.
%  ==========================================================================
if ~exist('harvey','var') || ~exist('maleL','var'), load(feasCheckpoint); end

kc4Rxns  = {'Colon_EX_ksii_core4(e)_[bc]','Liver_EX_ksii_core4(e)_[bc]','Kidney_EX_ksii_core4(e)_[bc]'};
kc4Stage = {'harvey (raw)', harvey; 'maleL (post-HMDB)', maleL};
fprintf('\nINV-3 — ksii_core4 uptake bound origin trace:\n');
for r = 1:numel(kc4Rxns)
    fprintf('  %-36s', kc4Rxns{r});
    for st = 1:size(kc4Stage,1)
        idx = strcmp(kc4Stage{st,2}.rxns, kc4Rxns{r});
        if any(idx)
            fprintf('  %-20s lb=%8.4f', kc4Stage{st,1}, kc4Stage{st,2}.lb(idx));
        end
    end
    fprintf('\n');
end
fprintf('  => Colon/Liver bounds from RAW model. Kidney relaxed by HMDB.\n');


%% ==========================================================================
%  INV-4.  Female choline shortfall: sIEC uptake cap too tight
%
%  ROADMAP:
%    INV-4a  Spleen pchol_hs: opening ub -> still stat=0, residual -1.2731.  NOT THE FIX.
%    INV-4b  Shadow prices: value 1.0 on sIEC Kennedy pathway
%            (chol->cholp->cdpchol->pchol_hs). Intestinal PC synthesis is the sink.
%    INV-4c  BBB cycling: capping BBB_CHOL[CSF]upt -1 to -500 leaves residual
%            invariant at -1.2731.  DECOUPLED — NOT THE FIX.
%    INV-4d  CAUSE: sIEC choline demand is 16.2 mmol/day (female) vs 14.0 (male).
%            sIEC_EX_chol(e)_[bc] was capped at -14.93, leaving a 1.28 shortfall.
%            Blood choline is abundant (~509 mmol/day released from organs).
%            Closing sIEC_EX_chol breaks even the male — it is the choline valve.
%    INV-4e  Binary search: minimum sIEC_EX_chol(e)_[bc] for female = -16.20 mmol/day.
%            Fix applied in OLD script: -20 (headroom). New approach: see INV-6.
%  ==========================================================================
if ~exist('femaleFastInfeasible','var'), load(feasCheckpoint); end

% INV-4a — spleen
if ~exist('inv4a','var')
    femBase = changeRxnBounds(femaleFastInfeasible, 'Muscle_sink_glu_L(c)', -280, 'l');
    s = optimizeWBModel(changeRxnBounds(femBase, 'Spleen_EX_pchol_hs(e)_[bc]', 1e6, 'u'));
    fprintf('\nINV-4a — spleen pchol_hs: stat=%d  (expected 0 — spleen is not the fix)\n', s.stat);
    inv4a = true;
end

% INV-4c — BBB cycling
if ~exist('inv4c','var')
    femBase = changeRxnBounds(femaleFastInfeasible, 'Muscle_sink_glu_L(c)', -280, 'l');
    fprintf('INV-4c — BBB_CHOL[CSF]upt cap sweep:\n');
    for cap = [-1 -5 -50 -500]
        s = optimizeWBModel(changeRxnBounds(femBase, 'BBB_CHOL[CSF]upt', cap, 'l'));
        fprintf('   lb=%6.0f: stat=%d  (all expected 0 — BBB is decoupled)\n', cap, s.stat);
    end
    inv4c = true;
end

% INV-4e — minimum sIEC choline bound
if ~exist('siecCholMin','var')
    femBase = changeRxnBounds(femaleFastInfeasible, 'Muscle_sink_glu_L(c)', -280, 'l');
    lo = -50; hi = -14.93;
    for it = 1:14
        mid = (lo + hi)/2;
        s = optimizeWBModel(changeRxnBounds(femBase, 'sIEC_EX_chol(e)_[bc]', mid, 'l'));
        if s.stat == 1, lo = mid; else, hi = mid; end
    end
    siecCholMin = lo;
end
fprintf('INV-4e — min sIEC_EX_chol(e)_[bc] for female feasibility: %.3f mmol/day\n', siecCholMin);
fprintf('         (= -14.93 cap + 1.27 shortfall = ~-16.20)\n');


%% ==========================================================================
%  INV-5.  noTrTr objective — a documented limitation
%
%  biomass_maintenance_noTrTr (no transcription/translation) is Thiele 2020's
%  intended fasting objective. Only 4 organs in Harvey have it (Platelet, RBC,
%  Testis, Prostate) and only 2 in Harvetta. Major organs (Muscle, Brain, Liver,
%  Kidney, Adipocytes) have no noTrTr variant.
%
%  Whole_body_objective_rxn is hardwired to biomass_maintenance_dummy_objective.
%  Switching available organs to noTrTr gives stat=0 — rewiring the objective
%  stoichiometry would be required.
%
%  OPEN QUESTION (for Ines): was this incomplete by design, or is it unfinished?
%  Could noTrTr variants be added to major organs and the objective rewired?
%  ==========================================================================
fprintf('\nINV-5 — noTrTr objective: only 4 (Harvey) / 2 (Harvetta) organs have it.\n');
fprintf('         Switching gives stat=0. Using biomass_maintenance as fasting objective.\n');
fprintf('         Ask Ines: incomplete implementation or deliberate biological choice?\n');


%% ==========================================================================
%  INV-6.  Bare minimum feasibility — what does the fasting model structurally need?
%
%  QUESTION: if storage sinks are unconstrained, what is the structural minimum
%  needed for feasibility (ions alone? choline alone? both?)?
%
%  METHOD: set fasting sink value to large negative number, add things one at a time.
%
%  RESULTS (both sexes):
%    Step 1 — bare fasting (-1e6 sinks, no diet):          male=0  female=0
%    Step 2 — + sIEC choline (-14.93), no ions:            male=0  female=0
%    Step 3 — + ions only, no choline:                     male=0  female=0
%    Step 4 — + ions + sIEC choline:                       male=1  female=1
%  => ions AND choline are both required; neither is sufficient alone.
%
%  THRESHOLD SEARCH (multiples of 10, exact values TBD):
%    Female: feasible at fasting value -20 (infeasible at -10)
%    Male:   feasible at fasting value -30 (infeasible at -20)
%  => male needs a higher sink release rate than female (larger muscle mass)
%
%  KEY INSIGHT: the GAG proxy, glutamate relaxation, and pancreas choline fixes
%  from the OLD script are all artefacts of the -10 sink cap. With unconstrained
%  (or less tight) sinks the model routes around those gaps internally.
%
%  IMPLICATION FOR fastingWBM.m: use sink values ~-20/-30 + ions + sIEC choline.
%  The -10 default is physiologically unjustified; the old fixes were compensating
%  for it. TODO: find exact threshold via binary search and justify with literature.
%  ==========================================================================
if ~exist('maleL','var') || ~exist('femaleL','var'), load(feasCheckpoint); end
if ~exist('fastingIons','var')
    fastingIons = {'Diet_EX_h2o[d]','Diet_EX_na1[d]','Diet_EX_k[d]','Diet_EX_pi[d]', ...
                   'Diet_EX_cl[d]','Diet_EX_so4[d]','Diet_EX_zn2[d]','Diet_EX_ca2[d]', ...
                   'Diet_EX_cu2[d]','Diet_EX_mn2[d]','Diet_EX_fe2[d]','Diet_EX_fe3[d]'};
end

% Step 1 — bare fasting
mD = setFeedingFastingConstraints(maleL,   'fasting', -1e6, 1e6); mD.osenseStr = 'max';
fD = setFeedingFastingConstraints(femaleL, 'fasting', -1e6, 1e6); fD.osenseStr = 'max';
s = optimizeWBModel(mD); fprintf('\nINV-6 step 1 — bare (-1e6, no diet):    male stat=%d\n', s.stat);
s = optimizeWBModel(fD); fprintf('INV-6 step 1 — bare (-1e6, no diet):  female stat=%d\n', s.stat);

% Step 2 — choline only
s = optimizeWBModel(changeRxnBounds(mD, 'sIEC_EX_chol(e)_[bc]', -14.93, 'l'));
fprintf('INV-6 step 2 — + choline only:          male stat=%d\n', s.stat);
s = optimizeWBModel(changeRxnBounds(fD, 'sIEC_EX_chol(e)_[bc]', -14.93, 'l'));
fprintf('INV-6 step 2 — + choline only:        female stat=%d\n', s.stat);

% Step 3 — ions only
s = optimizeWBModel(changeRxnBounds(mD, fastingIons, -1e6, 'l'));
fprintf('INV-6 step 3 — + ions only:             male stat=%d\n', s.stat);
s = optimizeWBModel(changeRxnBounds(fD, fastingIons, -1e6, 'l'));
fprintf('INV-6 step 3 — + ions only:           female stat=%d\n', s.stat);

% Step 4 — ions + choline
mD4 = changeRxnBounds(changeRxnBounds(mD, fastingIons, -1e6, 'l'), 'sIEC_EX_chol(e)_[bc]', -14.93, 'l');
fD4 = changeRxnBounds(changeRxnBounds(fD, fastingIons, -1e6, 'l'), 'sIEC_EX_chol(e)_[bc]', -14.93, 'l');
s = optimizeWBModel(mD4); fprintf('INV-6 step 4 — + ions + choline:        male stat=%d\n', s.stat);
s = optimizeWBModel(fD4); fprintf('INV-6 step 4 — + ions + choline:      female stat=%d\n', s.stat);

clearvars mD fD mD4 fD4 s
