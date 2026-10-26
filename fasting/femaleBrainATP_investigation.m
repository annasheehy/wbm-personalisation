%% ==========================================================================
%  femaleBrainATP_investigation.m
%
%  QUESTION.  Harvetta's fasting brain can supply only ~1514 mmol ATP/day at
%  obj=1 — below the 3500 resting floor — so the female fasting BMF is infeasible.
%  Harvey's brain reaches ~3636. Both hearts have ATP to spare (F 8724 / M 9457
%  vs the 1000 floor), so this is NOT a whole-body energy shortage: the brain
%  cannot oxidise fatty acids (blood-brain barrier) and runs on glucose + ketones,
%  so the female brain must be short of those specific fuels.
%
%  APPROACH (most decisive first).
%    1  Re-establish the brain-ATP ceilings for both sexes; capture the female
%       solution (fluxes + reduced costs).
%    2  CAUSAL knob test — give the FEMALE model the MALE fasting parameters
%       (glycogen -969, storage-AA sinks -30) and re-max the brain. Recovery =>
%       the bottleneck is our fasting knobs (female-specific choices, fixable);
%       no recovery => female physiology (HMDB bounds / organ scaling), deeper.
%    3  Reduced-cost bottleneck — reactions AT a bound that limit the female brain
%       max: the LP's own answer to "what is binding", no guessing.
%    4  Fuel divergence — brain glucose / ketone uptake, male vs female.
%
%  Assumes maleFast / femaleFast / maleL / femaleL in the workspace (run
%  fastingWBM.m Parts 0-3 first). Read-only except relaxation tests (local copies).
%  Reaction names found below by pattern match are flagged 'confirm' — NEVER
%  assumed correct; search model.rxns / ask before concluding one is absent.
%  ==========================================================================
%
%  --------------------------------------------------------------------------
%  FINDINGS LOG (9 Sep 2026) — what changes Harvetta's fasting brain ATP ceiling
%  Baseline (at AA -30): female brain max 3320 vs floor 3500 (SHORT 180); male 3636.
%
%  >>> UPDATE (9 Sep 2026): RESOLVED. The female fasting brain clears the 3500
%      floor once storage-AA sinks reach -40. The ceiling SCALES with AA magnitude
%      (sweep 2b): -20->1514, -30->3320, -40->4587 [CLEARS], -50->5474, -100->7420.
%      So the limiter is STORAGE-AA SUBSTRATE, not oxidation/delivery: at -50 the
%      female brain oxidises normally (O2 -3234, lactate export ~3) vs the -30
%      artefact (lactate export 307) that I misread as an oxidation limit below.
%      CHOSEN VALUE: -40 (minimal feasible; male already clears at -30). Floor-
%      scaling (section 8) NOT needed. OPEN ITEM: justify the AA-release magnitude
%      from literature. Numbers below are the -30 trail, kept for the record.
%
%  LEVERS THAT RAISE IT
%    * Storage-AA sinks -20 -> -30 : 1514 -> 3320  (DOMINANT; set by
%      setFeedingFastingConstraints). PREDICTED to plateau near 3320 beyond -30
%      (brain becomes oxidation- not substrate-limited) — tested by the sweep 2b.
%  LEVERS THAT DON'T
%    * Hepatic glycogen -754 -> -969 : NO change (1514->1514 at AA-20;
%      3320->3320 at AA-30). Glycogen release is not limiting the brain.
%    * Swap ALL 4027 female Brain_/BBB_ bounds -> male values : 3320 -> 3327
%      (+7, negligible) => the limit is NOT brain-internal or BBB reactions.
%    * Bounds-conflict ruled out: Brain_DM_atp_c_ is a SINGLE rxn, lb 0 / ub 1e6
%      in both sexes; 3320 is a genuine NETWORK maximum, not a ub cap.
%
%  MECHANISM (sections 4/5, at brain-max)
%    * Glucose capped at -100 in BOTH sexes (Brain_EX_glc_D(e)_[csf]).
%    * Female brain EXPORTS lactate +307 vs male IMPORTS -76 -> female ferments
%      the carbon overflow instead of oxidising it (oxidation-limited).
%    * Brain O2 uptake female -2066 vs male -2201; ATP synthase Brain_ATPS4m
%      female 8619 vs male 9376 (~8% less) -> tracks the ~9% lower ceiling.
%
%  ROOT (sections 6/7)
%    * Every BBB uptake bound is female ~0.927 x male (uniform organ-weight /
%      blood-flow scaling): O2, GLC_D, LAC_L, NA1 all 0.927. They have headroom
%      (not binding), so the limit is UPSTREAM of the brain: systemic blood
%      O2/fuel delivery = genuine Harvetta physiology (smaller organs & flows).
%
%  PROPOSED RESOLUTION (to test — section 8, TODO)
%    * The 3500 brain (and 1000 heart) ATP floor is FIXED/sex-invariant (Thiele
%      hardcoded = 3500 for both), while everything else about Harvetta is ~91%
%      scaled. Scale the electrophysiology floor by brain mass:
%      female floor ~ 3500 x (Harvetta brain wt / reference) ~ 3200 < 3320 =>
%      feasible. Fits the personalisation theme. Weights live in
%      IndividualParameters.OrgansWeights / calcOrganFract (OrganWeight[Fract]).
%
%  SECONDARY (section 3 reduced costs, for later)
%    * Glycan/GAG reactions dominate the shadow prices: Brain_HSAT2ly (164),
%      NACHEX/NS26Tg/GCNTg/CORE7GTg/A4GNT2g, Diet_EX_ksi (cf INV-3 missing
%      cartilage), and cardiolipin synthase CLS_hs (Muscle/Adipo/Liver ~93).
%      Competing biomass-maintenance demands, not the primary brain-ATP limit.
%  --------------------------------------------------------------------------
% Run standalone: prefer workspace vars (if you just ran fastingWBM.m), else load
% the Part-3 checkpoint directly — so you NEVER need to re-run the full pipeline.
if ~all([exist('maleFast','var') exist('femaleFast','var') exist('femaleL','var')])
    P = persWBMpaths();
    ck = fullfile(P.results, 'fasting', 'ckpt_feasibleModels.mat');
    assert(isfile(ck), 'Checkpoint not found (%s) — run fastingWBM.m Parts 0-3 first.', ck);
    load(ck, 'maleL', 'femaleL', 'maleFed', 'femaleFed', 'maleFast', 'femaleFast');
    fprintf('loaded fasting models from checkpoint (no pipeline re-run): %s\n', ck);
end

brainFloor  = 3500;
fastingIons = {'Diet_EX_h2o[d]','Diet_EX_na1[d]','Diet_EX_k[d]','Diet_EX_pi[d]', ...
               'Diet_EX_cl[d]','Diet_EX_so4[d]','Diet_EX_zn2[d]','Diet_EX_ca2[d]', ...
               'Diet_EX_cu2[d]','Diet_EX_mn2[d]','Diet_EX_fe2[d]','Diet_EX_fe3[d]'};

%% 0.1  Baseline brain-ATP ceilings + capture the female solution --------------
% fix obj=1, c(:)=0, maximise total Brain_DM_atp demand; plain LP (no norm).
capture = struct();
sexes = {'male', maleFast; 'female', femaleFast};
for q = 1:2
    sxn = sexes{q,1};  m = sexes{q,2};
    ob = strcmp(m.rxns, 'Whole_body_objective_rxn');  m.lb(ob) = 1;  m.ub(ob) = 1;
    bi = contains(m.rxns, 'Brain_DM_atp');
    m.c(:) = 0;  m.c(bi) = 1;  m.osenseStr = 'max';
    s = optimizeWBModel(m);
    capture.(sxn).model = m;  capture.(sxn).sol = s;  capture.(sxn).brainIdx = bi;
    fprintf('0.1 — %-6s brain ceiling %7.1f / floor %d  (stat %d)\n', sxn, s.f, brainFloor, s.stat);
end


%% 2  CAUSAL knob test — female physiology with male fasting knobs ------------
% Rebuild the female fasting model swapping in male values for the two knobs we
% chose per-sex (storage-AA sink magnitude, hepatic glycogen release), then
% re-max the brain. This isolates "our knobs" from "female physiology".
fprintf('\n2 — knob test: female model, brain ceiling after swapping in male knobs\n');

knob = {
   'female baseline (AA -20, glyc -754)',  -20, -754
   'female + male glycogen (-969)',        -20, -969
   'female + male AA sinks (-30)',         -30, -754
   'female + both male knobs',             -30, -969
};
for t = 1:size(knob,1)
    aaSink = knob{t,2};  glyc = knob{t,3};
    [~] = evalc('fm = setFeedingFastingConstraints(femaleL, ''fasting'', aaSink, abs(aaSink));');  % silence storageAA print + asn warning
    fm = changeRxnBounds(fm, fastingIons, -1e6, 'l');
    fm = changeRxnBounds(fm, 'sIEC_EX_chol(e)_[bc]', -6, 'l');
    fm = changeRxnBounds(fm, 'Liver_sink_glygn2(c)', glyc, 'l');
    fm = changeRxnBounds(fm, 'Liver_sink_glygn2(c)',    0, 'u');
    ob = strcmp(fm.rxns, 'Whole_body_objective_rxn');  fm.lb(ob) = 1;  fm.ub(ob) = 1;
    bi = contains(fm.rxns, 'Brain_DM_atp');
    fm.c(:) = 0;  fm.c(bi) = 1;  fm.osenseStr = 'max';
    s = optimizeWBModel(fm);
    flag = 'still SHORT';  if s.f >= brainFloor, flag = 'RECOVERS >= floor'; end
    fprintf('   %-38s brain %7.1f  (stat %d)  %s\n', knob{t,1}, s.f, s.stat, flag);
end


%% 2b  AA-sink sweep — does raising storage-AA beyond -30 lift the female brain?
% -30 was the plateau candidate in the knob test. The oxidation-limit mechanism
% (section 4) predicts more AA substrate will NOT help past ~-30 — test it, since
% AA sinks can be pushed to -100 as a generous fasting placeholder. Glycogen held
% at the female value (-754); only the storage-AA magnitude varies.
fprintf('\n2b — female brain ceiling vs storage-AA sink magnitude (glyc -754):\n');
for aa = [-20 -30 -40 -60 -100]
    [~] = evalc('fm = setFeedingFastingConstraints(femaleL, ''fasting'', aa, abs(aa));');
    fm = changeRxnBounds(fm, fastingIons, -1e6, 'l');
    fm = changeRxnBounds(fm, 'sIEC_EX_chol(e)_[bc]', -6, 'l');
    fm = changeRxnBounds(fm, 'Liver_sink_glygn2(c)', -754, 'l');
    fm = changeRxnBounds(fm, 'Liver_sink_glygn2(c)',    0, 'u');
    ob = strcmp(fm.rxns, 'Whole_body_objective_rxn');  fm.lb(ob) = 1;  fm.ub(ob) = 1;
    bi = contains(fm.rxns, 'Brain_DM_atp');  fm.c(:) = 0;  fm.c(bi) = 1;  fm.osenseStr = 'max';
    s = optimizeWBModel(fm);
    flag = 'SHORT';  if s.f >= brainFloor, flag = 'CLEARS floor'; end
    fprintf('   AA %4d : brain %7.1f  (stat %d)  %s\n', aa, s.f, s.stat, flag);
end


%% 3  Reduced-cost bottleneck at the female brain-max ------------------------
% Reactions sitting AT a bound with nonzero reduced cost are what limits the
% female brain objective. This is the LP's own answer — no name guessing.
fprintf('\n3 — female brain-max: top binding reactions (at a bound, |reduced cost|>0)\n');
sF = capture.female.sol;  mF = capture.female.model;
if isfield(sF, 'w') && ~isempty(sF.w) && ~isempty(sF.v)
    tol     = 1e-6;
    atLB    = abs(sF.v - mF.lb) < tol;
    atUB    = abs(sF.v - mF.ub) < tol;
    binding = (atLB | atUB) & abs(sF.w) > tol;
    idx     = find(binding);
    [~, ord] = sort(abs(sF.w(idx)), 'descend');
    idx = idx(ord);
    fprintf('   %-34s %10s %10s %10s %10s\n', 'reaction', 'v', 'lb', 'ub', 'redcost');
    for k = 1:min(25, numel(idx))
        j = idx(k);
        fprintf('   %-34s %10.2f %10.1f %10.1f %10.3f\n', ...
            mF.rxns{j}, full(sF.v(j)), mF.lb(j), mF.ub(j), sF.w(j));
    end
    if isempty(idx), fprintf('   (no binding reactions with nonzero reduced cost found)\n'); end
else
    fprintf('   reduced costs (solution.w) not returned — rerun, or inspect duals (solution.y).\n');
end


%% 4  Fuel divergence — brain glucose / ketone uptake, male vs female --------
% Candidate brain fuel reactions by pattern (CONFIRM these names — do not assume
% a fuel is absent from a non-match). Male/female reaction indices differ
% (sex-specific organs), so female reactions are matched to male BY NAME.
fprintf('\n4 — brain fuel candidates (names to CONFIRM) and flux at brain-max:\n');
fuelPat = {'glucose','glc_D'; 'ketone-bhb','bhb'; 'ketone-acac','acac'; ...
           'lactate','lac_L'; 'oxygen','o2'};
vM = capture.male.sol.v;   vF = capture.female.sol.v;
for f = 1:size(fuelPat,1)
    tokRe = ['_' fuelPat{f,2} '[\(\[]'];              % anchored: id then ( or [ (excludes co2/h2o2)
    ids   = find(contains(maleFast.rxns, 'Brain_') & ...
                 ~cellfun(@isempty, regexp(maleFast.rxns, tokRe, 'once')));
    mvv  = nan(numel(ids),1);  fvv = nan(numel(ids),1);  names = cell(numel(ids),1);  cnt = 0;
    for k = 1:numel(ids)
        rn = maleFast.rxns{ids(k)};
        jF = find(strcmp(femaleFast.rxns, rn), 1);
        fv = NaN;  if ~isempty(jF), fv = full(vF(jF)); end
        mv = full(vM(ids(k)));
        if abs(mv) > 1e-3 || (~isnan(fv) && abs(fv) > 1e-3)
            cnt = cnt + 1;  mvv(cnt) = mv;  fvv(cnt) = fv;  names{cnt} = rn;
        end
    end
    mvv = mvv(1:cnt);  fvv = fvv(1:cnt);  names = names(1:cnt);
    fprintf('   [%s] %d brain "%s" rxns carry flux:\n', fuelPat{f,1}, cnt, fuelPat{f,2});
    [~, o] = sort(abs(fvv), 'descend', 'MissingPlacement', 'last');
    for k = 1:min(6, cnt)
        fprintf('       male %10.2f   female %10.2f   %s\n', mvv(o(k)), fvv(o(k)), names{o(k)});
    end
end


%% 5  Brain ATP-producing reactions, male vs female (indicative) -------------
% Net ATP production per BRAIN reaction at brain-max, computed per sex with that
% sex's OWN indices (male/female reaction sets differ). Alternative optima make
% raw flux comparison indicative, not definitive — read alongside 2 and 3.
fprintf('\n5 — top brain ATP-producing reactions at brain-max (net atp made):\n');
atpMetsM = find(contains(maleFast.mets,   '_atp['));   brainRxM = find(contains(maleFast.rxns,   'Brain_'));
atpMetsF = find(contains(femaleFast.mets, '_atp['));   brainRxF = find(contains(femaleFast.rxns, 'Brain_'));
netM = full(sum(maleFast.S(atpMetsM, brainRxM)   .* vM(brainRxM)', 1));   % >0 = net producer
netF = full(sum(femaleFast.S(atpMetsF, brainRxF) .* vF(brainRxF)', 1));
[~, oM] = sort(netM, 'descend');   [~, oF] = sort(netF, 'descend');
fprintf('   -- male --\n');
for k = 1:8, j = brainRxM(oM(k)); fprintf('   %10.1f  %s\n', netM(oM(k)), maleFast.rxns{j}); end
fprintf('   -- female --\n');
for k = 1:8, j = brainRxF(oF(k)); fprintf('   %10.1f  %s\n', netF(oF(k)), femaleFast.rxns{j}); end


%% 6  Which Brain/BBB reactions are more tightly bounded in female than male? -
% The user's original question, now targeted. Name-matched (male/female indices
% differ). "Tighter" = female allows less flux in the direction it would need
% (higher lb => less uptake; lower ub => less output). Needs section 0.1 run.
fprintf('\n6 — Brain/BBB reactions more tightly bounded in female than male:\n');
brx  = find(contains(maleFast.rxns, 'Brain_') | contains(maleFast.rxns, 'BBB'));
gap  = nan(numel(brx),1);  info = cell(numel(brx),1);  cnt = 0;
for k = 1:numel(brx)
    jF = find(strcmp(femaleFast.rxns, maleFast.rxns{brx(k)}), 1);
    if isempty(jF), continue; end
    dLB = femaleFast.lb(jF) - maleFast.lb(brx(k));   % >0: female lb higher -> less uptake capacity
    dUB = maleFast.ub(brx(k)) - femaleFast.ub(jF);   % >0: female ub lower  -> less output capacity
    g   = max(dLB, dUB);
    if g > 1e-6
        cnt = cnt + 1;  gap(cnt) = g;
        info{cnt} = sprintf('%-30s  F[%g, %g]   M[%g, %g]', maleFast.rxns{brx(k)}, ...
            femaleFast.lb(jF), femaleFast.ub(jF), maleFast.lb(brx(k)), maleFast.ub(brx(k)));
    end
end
gap = gap(1:cnt);  info = info(1:cnt);
[~, o] = sort(gap, 'descend');
for k = 1:min(25, cnt), fprintf('   %s\n', info{o(k)}); end
if cnt == 0, fprintf('   (no Brain/BBB reaction is tighter in female)\n'); end


%% 7  Causal swap — give the female brain MALE Brain/BBB bounds, then re-max ---
% If the female ceiling rises toward the male 3636 when its Brain_/BBB_ bounds are
% set to male values, the limit is brain-internal/BBB (trace via section 6). If it
% stays ~3320, the limit is UPSTREAM (blood fuel/O2 delivery) -> HMDB blood bounds
% or organ-weight scaling, not the brain itself. Needs section 0.1 run.
fprintf('\n7 — swap female Brain/BBB bounds -> male values, then re-max brain:\n');
fSwap = femaleFast;  nSw = 0;
bswap = find(contains(maleFast.rxns, 'Brain_') | contains(maleFast.rxns, 'BBB'));
for k = 1:numel(bswap)
    jF = find(strcmp(fSwap.rxns, maleFast.rxns{bswap(k)}), 1);
    if isempty(jF), continue; end
    fSwap.lb(jF) = maleFast.lb(bswap(k));  fSwap.ub(jF) = maleFast.ub(bswap(k));  nSw = nSw + 1;
end
ob = strcmp(fSwap.rxns, 'Whole_body_objective_rxn');  fSwap.lb(ob) = 1;  fSwap.ub(ob) = 1;
bi = contains(fSwap.rxns, 'Brain_DM_atp');  fSwap.c(:) = 0;  fSwap.c(bi) = 1;  fSwap.osenseStr = 'max';
sSw = optimizeWBModel(fSwap);
fprintf('   swapped %d Brain/BBB bounds; female brain now %7.1f  (was %.1f, male %.1f, floor %d)\n', ...
    nSw, sSw.f, capture.female.sol.f, capture.male.sol.f, brainFloor);
if sSw.f >= brainFloor
    fprintf('   -> RECOVERS: limit is brain-internal / BBB bounds (trace via section 6).\n');
else
    fprintf('   -> still short: limit is UPSTREAM (blood fuel / O2 delivery), not the brain itself.\n');
end
