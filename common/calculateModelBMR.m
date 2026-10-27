function bmr = calculateModelBMR(model, varargin)
% calculateModelBMR  Basal metabolic rate (BMR) of a whole-body model (WBM).
%
%   bmr = calculateModelBMR(model)
%   bmr = calculateModelBMR(model, 'Name', Value, ...)
%
% Computes the model's BMR (kcal/day) by the Thiele 2020 basal-metabolic-flux
% method (papers/2020_WholeBodyModelling/MethodSection3.mlx §3.3.2.1):
%   1. BMR conditions: resting brain and heart ATP demands, whole-body maintenance
%      fixed to 1 (setBMRConditions);
%   2. the basal flux distribution: the flux vector with the smallest Euclidean norm
%      (QP, min ||v||_2);
%   3. ATP consumption summed over all organs and compartments and converted to kcal
%      (bmrFromFlux; 64 kJ/mol ATP).
%
% The plain QP does not converge on the ~80k-reaction WBM (Gurobi stat -1) because of
% badly scaled coupling constraints. It converges after liftCouplingConstraints (Sun,
% Fleming, Thiele 2013) and a barrier-only solve with crossover off.
%
% Name-Value options (defaults):
%   'brainATPfloor'  3500   lb on Brain_DM_atp (mmol/day)
%   'heartATPfloor'  1000   lb on Heart_DM_atp (mmol/day)
%   'kJperATP'       64     energy per mol ATP hydrolysed (Wackerhage 1998)
%   'minNorm'        1e-6   weight of the 2-norm (paper value). With no linear objective the
%                           solution does not depend on it; it is increased only if a solve
%                           does not converge
%   'timeLimit'      1800   solver time limit (s)
%   'bigCoeff'       1e4    BIG threshold for liftCouplingConstraints (skipped if lifted)
%
% OUTPUT (struct): kcal  BMR (kcal/day), NaN if no solution
%                  atp   ATP consumed (mmol/day), NaN if no solution
%                  stat  solver status (1 optimal; 3 suboptimal; otherwise no solution)
%                  v     basal flux distribution, original reactions only (empty if none)
%
% Anna Sheehy 2026.

p = inputParser;
addParameter(p, 'brainATPfloor', 3500, @isnumeric);
addParameter(p, 'heartATPfloor', 1000, @isnumeric);
addParameter(p, 'kJperATP',      64,   @isnumeric);
addParameter(p, 'minNorm',       1e-6, @isnumeric);
addParameter(p, 'timeLimit',     1800, @isnumeric);
addParameter(p, 'bigCoeff',      1e4,  @isnumeric);
parse(p, varargin{:});
opt = p.Results;

m = setBMRConditions(model, opt.brainATPfloor, opt.heartATPfloor);
m.c(:) = 0;  m.osenseStr = 'min';
if ~isModelLifted(m) && isfield(m, 'C') && ~isempty(m.C)
    m = liftCouplingConstraints(m, opt.bigCoeff, 0);
end
param = struct('minNorm', opt.minNorm, 'secondsTimeLimit', opt.timeLimit);
if strcmpi(getCobraSolver('QP'), 'gurobi')
    param.Method = 2;  param.Crossover = 0;  param.NumericFocus = 3;  param.BarHomogeneous = 1;
    param.TimeLimit = opt.timeLimit;
end
sol = struct('stat', -1, 'v', []);
for lam = opt.minNorm * [1 1e2 1e4]                 % larger weight only if a solve fails
    param.minNorm = lam;
    sol = optimizeWBModel(m, param);
    if isfield(sol, 'v') && ~isempty(sol.v) && ismember(sol.stat, [1 3]), break, end
end

bmr = struct('kcal', NaN, 'atp', NaN, 'stat', sol.stat, 'v', []);
if ismember(sol.stat, [1 3]) && isfield(sol, 'v') && ~isempty(sol.v)
    bmr.v = sol.v(1:numel(model.rxns));             % lifted solves append auxiliary variables
    [bmr.kcal, bmr.atp] = bmrFromFlux(model, bmr.v, opt.kJperATP);
end
end
