function model = setFastingConstraints(model, sex, varargin)
% setFastingConstraints  Apply a complete, feasible fasting constraint set to a WBM.
%
%   model = setFastingConstraints(model, sex)
%   model = setFastingConstraints(model, sex, 'Name', Value, ...)
%
% Wraps setFeedingFastingConstraints(model,'fasting',...) (the storage
% fill<->drain switch) and layers on the additional constraints needed for a
% FEASIBLE fasting whole-body model:
%   1. storage release enabled (organ sinks lb<0, ub=0)  [setFeedingFastingConstraints]
%      and dietary intake closed;
%   2. dietary ions capped at EU-average-diet provision (all except phosphate);
%   2b. endogenous phosphate (bone Pi — the model has no bone organ);
%   3. hepatic glycogen release enabled (sex-specific, no synthesis);
%   4. endogenous intestinal choline uptake from blood (sex-specific).
%
% This turns "one function call + several scattered manual steps" into one
% reproducible call. It does not modify setFeedingFastingConstraints.
%
% INPUT
%   model   whole-body model (e.g. Harvey / Harvetta)
%   sex     'male' | 'female'  (selects the physiological defaults below)
%
% Name-Value options (defaults):
%   'aaRelease'  -40         storage-AA / organ-store release lb (the fastingValue
%                            passed to setFeedingFastingConstraints). PROVISIONAL.
%   'glycogen'   sex         Liver_sink_glygn2(c) release lb; default male -969 /
%                            female -754 (Jensen 2011, Olson 1984). ub set to 0.
%   'choline'    sex         sIEC_EX_chol(e)_[bc] uptake lb from blood; default
%                            male -2 / female -6.
%   'ionCaps'    {11x2}      Diet_EX_<ion>[d] uptake caps {rxn, lb}; lb = -EU-average
%                            provision (mmol/day; from EUAverageDietNew). All 11
%                            non-phosphate ions; Part 5 confirmed none bind under
%                            the resting brain floor.
%   'phosphate'  -1e6        Diet_EX_pi[d] lb — endogenous phosphate release (bone/
%                            tissue) during fasting; NOT dietary. REQUIRED to meet the
%                            resting brain ATP floor: the model has no bone organ, so
%                            capping Pi at dietary ~14 mmol/day starves the brain
%                            (max Brain_DM_atp ~2600 < 3500). [TODO: justify the rate]
%
% OUTPUT
%   model    with fasting constraints applied; model.SetupInfo.FeedingStatus='fasting'
%            and model.SetupInfo.FastingParams recording the values used.
%
% NOTE: the release magnitudes are PROVISIONAL placeholders chosen for feasibility
% (aaRelease -40 is the minimum at which the female brain meets the 3500 resting-ATP
% floor: female brain max -20 -> 1514, -30 -> 3320, -40 -> 4587). They still need
% literature justification (cf Felig 1970; see femaleBrainATP_investigation.m).
%
% Anna Sheehy 2026.  Companion to setFeedingFastingConstraints (Ines Thiele).

sex = validatestring(lower(sex), {'male','female'});

p = inputParser;
addParameter(p, 'aaRelease', -40,  @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'glycogen',  NaN,  @(x) isnumeric(x) && isscalar(x));   % NaN -> sex default
addParameter(p, 'choline',   NaN,  @(x) isnumeric(x) && isscalar(x));   % NaN -> sex default
addParameter(p, 'ionCaps', { ...     % lb = -EU-average-diet provision (mmol/day), from EUAverageDietNew
    'Diet_EX_h2o[d]', -165892.3425;  'Diet_EX_na1[d]', -172.4463;  'Diet_EX_k[d]',  -83.4338; ...
    'Diet_EX_cl[d]',    -187.7787;   'Diet_EX_so4[d]',        0;    'Diet_EX_zn2[d]',  -0.1510; ...
    'Diet_EX_ca2[d]',    -18.5491;   'Diet_EX_cu2[d]',   -0.0357;   'Diet_EX_mn2[d]',  -0.0747; ...
    'Diet_EX_fe2[d]',     -0.1092;   'Diet_EX_fe3[d]',   -0.1092}, @(x) iscell(x) && size(x,2) == 2);
addParameter(p, 'phosphate', -1e6, @(x) isnumeric(x) && isscalar(x));   % endogenous Pi (bone)
parse(p, varargin{:});
opt = p.Results;

% sex-specific defaults, only where the caller did not override them
if isnan(opt.glycogen), opt.glycogen = pickSex(sex, -969, -754); end
if isnan(opt.choline),  opt.choline  = pickSex(sex,   -2,   -6); end

% 1. storage flip: enable release from organ stores, close dietary intake.
%    (In 'fasting' mode setFeedingFastingConstraints ignores its 4th arg, so the
%    storage upper bound is fixed at 0 there; abs(aaRelease) is a harmless filler.)
model = setFeedingFastingConstraints(model, 'fasting', opt.aaRelease, abs(opt.aaRelease));

% 2. dietary ions capped at EU-average-diet provision (from EUAverageDietNew). The
%    Part-5 investigation confirmed none of these bind under the resting brain floor.
%    PHOSPHATE is the exception -> endogenous (handled in 2b).
ionRxns = opt.ionCaps(:,1);
present = ismember(ionRxns, model.rxns);
if any(~present)
    warning('setFastingConstraints:missingIons', '%d ion exchange(s) not in model, skipped: %s', ...
        sum(~present), strjoin(ionRxns(~present), ', '));
end
for i = find(present(:)')
    model = changeRxnBounds(model, opt.ionCaps{i,1}, opt.ionCaps{i,2}, 'l');
end

% 2b. endogenous phosphate. Bone stores ~85% of body phosphate and mobilises it
%     during fasting, but Harvey/Harvetta have NO bone organ (Thiele 2020), so Pi
%     must enter through Diet_EX_pi[d]. This is NOT dietary phosphate — it stands in
%     for bone/tissue Pi release, and is REQUIRED to meet the resting brain ATP floor
%     (capping Pi at dietary ~14 mmol/day starves it: max Brain_DM_atp ~2600 < 3500).
model = setRequired(model, 'Diet_EX_pi[d]', opt.phosphate, []);

% 3. hepatic glycogen release, no synthesis (lb<0, ub=0)
model = setRequired(model, 'Liver_sink_glygn2(c)', opt.glycogen, 0);

% 4. endogenous intestinal choline uptake from blood (lb only)
model = setRequired(model, 'sIEC_EX_chol(e)_[bc]', opt.choline, []);

% provenance
model.SetupInfo.FeedingStatus = 'fasting';
model.SetupInfo.FastingParams = opt;
end


function v = pickSex(sex, maleVal, femaleVal)
if strcmp(sex, 'male'), v = maleVal; else, v = femaleVal; end
end


function model = setRequired(model, rxn, lb, ub)
% Set bounds on a REQUIRED named reaction; warn clearly if it is absent, rather
% than failing silently (the fasting setup depends on it — flag a wrong/renamed id).
if ~any(strcmp(model.rxns, rxn))
    warning('setFastingConstraints:rxnMissing', ...
        'Reaction "%s" not found — fasting constraint NOT applied. Check the name against model.rxns.', rxn);
    return
end
if ~isempty(lb), model = changeRxnBounds(model, rxn, lb, 'l'); end
if ~isempty(ub), model = changeRxnBounds(model, rxn, ub, 'u'); end
end
