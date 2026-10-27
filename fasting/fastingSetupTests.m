%% fastingSetupTests.m
% Thesis Ch5, fasted state: the two tests that set the fasted constraints of
% setFastingConstraints. For each, the largest brain and heart ATP supply of the fasted
% models (LP, whole-body maintenance fixed to 1) is compared with the resting demands
% used for BMR (brain 3500, heart 1000 mmol/day).
%   1. release limit of each amino-acid store: 20, 30, 40 mmol/day (aaRelease)
%   2. phosphate uptake capped at its dietary provision (EUAverageDietNew) or left open
% Uses the lifted base models of fastingWBM.m (ckpt_liftedBase.mat); saves
% results/fasting/fastingSetupTests.csv. NaN = infeasible
% (e.g. Harvey at 20 mmol/day). Results: Harvetta brain 1514 / 3320 / 4587 at 20 / 30 / 40,
% Harvey 3636 at 30; Harvetta with phosphate capped at 14.2 mmol/day 2464, open 4587.
%
% Anna Sheehy 2026.

P = persWBMpaths();
resPath = fullfile(P.results, 'fasting');
global CBT_LP_SOLVER
if isempty(CBT_LP_SOLVER), initCobraToolbox(false); end
changeCobraSolver('gurobi', 'LP', 0);
Lf = load(fullfile(resPath, 'ckpt_liftedBase.mat'), 'maleL', 'femaleL');
base = struct('male', Lf.maleL, 'female', Lf.femaleL);
EUAverageDietNew;                                                    % defines Diet
dietPi = -str2double(Diet{strcmp(Diet(:, 1), 'Diet_EX_pi[d]'), 2});  % mmol/day

out = cell(0, 5);                                                    % sex, aaRelease, phosphate limit, brain, heart
for s = {'male', 'female'}
    for aa = [-20 -30 -40]
        [brain, heart] = maxOrganATP(setFastingConstraints(base.(s{1}), s{1}, 'aaRelease', aa));
        fprintf('%-6s amino-acid release %d: brain ATP max %.0f, heart ATP max %.0f mmol/day\n', ...
            s{1}, -aa, brain, heart);
        out(end+1, :) = {s{1}, -aa, Inf, brain, heart}; %#ok<SAGROW>
    end
end
for piLimit = [dietPi -1e6]
    [brain, heart] = maxOrganATP(setFastingConstraints(base.female, 'female', 'phosphate', piLimit));
    fprintf('female phosphate uptake limit %g: brain ATP max %.0f, heart ATP max %.0f mmol/day\n', ...
        -piLimit, brain, heart);
    out(end+1, :) = {'female', 40, -piLimit, brain, heart}; %#ok<SAGROW>
end
T = cell2table(out, 'VariableNames', {'sex', 'aaRelease', 'phosphateLimit', 'brainATPmax', 'heartATPmax'});
T.phosphateLimit(T.phosphateLimit >= 1e6) = Inf;                     % open
writetable(T, fullfile(resPath, 'fastingSetupTests.csv'));

function [brain, heart] = maxOrganATP(m)
% largest brain and heart ATP hydrolysis (mmol/day) with whole-body maintenance fixed to 1;
% NaN if the model is infeasible
m = setBMRConditions(m, [], []);
f = nan(1, 2);  org = {'Brain_DM_atp', 'Heart_DM_atp'};
for k = 1:2
    mk = m;  mk.c(:) = 0;  mk.c(contains(mk.rxns, org{k})) = 1;  mk.osenseStr = 'max';
    sol = optimizeWBModel(mk);
    if sol.stat == 1, f(k) = sol.f; end
end
brain = f(1);  heart = f(2);
end
