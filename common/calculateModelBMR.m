function bmr = calculateModelBMR(model, varargin)
% calculateModelBMR  Basal metabolic rate (BMR) of a whole-body model (WBM).
%
%   bmr = calculateModelBMR(model)
%   bmr = calculateModelBMR(model, 'Name', Value, ...)
%
% Computes the model's BMR (kcal/day) by the Thiele 2020 basal-metabolic-flux
% method (papers/2020_WholeBodyModelling/MethodSection3.m §3.3.2.1):
%   1. impose resting energy floors (brain electrophysiology, cardiac pumping);
%   2. find the parsimonious resting flux distribution (minimise a flux norm,
%      which strips thermodynamically infeasible futile cycles);
%   3. sum ATP hydrolysis across ALL compartments/organs;
%   4. convert: kcal = ATP(mmol/day) * kJperATP / 1000 * 0.239006  (1 kJ = 0.239006 kcal).
%
% Two norms are reported:
%   L2 : minimise ||v||_2 (QP) — the norm used in MethodSection3. The plain QP is
%        numerically unstable on the ~80k-reaction WBM (Gurobi -> stat=-1; CPLEX ->
%        stat=3) because badly-scaled COUPLING constraints wreck the barrier. This is
%        fixed by liftCouplingConstraints (Sun, Fleming, Thiele 2013: reformulate the
%        coupling into a better-scaled, higher-dimensional problem) followed by a
%        barrier-only solve with crossover OFF (stat=3 is usually crossover failing to
%        recover a vertex from an already-optimal barrier point). Together these give
%        stat=1 at the paper's minNorm=1e-6, for both sexes.
%   L1 : minimise sum|v|  (LP) via optimizeCbModel(m, 'min', 'one'). L1 CANNOT be run
%        through optimizeWBModel: it forwards minNorm='one' into solveCobraLP, whose
%        `minNorm ~= 0 &&` test errors on a non-scalar minNorm. optimizeCbModel handles
%        'one' internally, and still honours the coupling constraints (which
%        buildOptProblemFromModel appends as A=[S;C], b=[b;d]).
%
% L2 is the validated headline (reproduces the Thiele archetypes); L1 is a sparse
% cross-check that under-reads total ATP turnover (~half of L2).
%
% Name-Value options (defaults):
%   'brainATPfloor'  3500        lb on reactions matching 'Brain_DM_atp' (mmol/day)
%   'heartATPfloor'  1000        lb on reactions matching 'Heart_DM_atp' (mmol/day)
%   'fixObjective'   true        §3.3 fixes Whole_body_objective_rxn = 1 before the
%                                BMF solve (MethodSection3 lines 769-770:
%                                changeRxnBounds(...,1,'b')); true matches the paper.
%   'kJperATP'       64          energy per mol ATP hydrolysed (Wackerhage 1998)
%   'minNorm'        1e-6        2-norm weight for the L2 QP (paper value). With c=0 the
%                                solution is invariant to it; it only affects conditioning.
%                                Escalated automatically if a solve does not converge.
%   'timeLimit'      1800        secondsTimeLimit for the L2 optimizeWBModel solve
%   'lift'           true        reformulate coupling constraints (liftCouplingConstraints)
%                                before the L2 QP — the numerical fix for WBM scale
%   'bigCoeff'       1e4         BIG threshold passed to liftCouplingConstraints (1e3-1e4)
%   'norms'          {'L1','L2'} which norm(s) to compute
%
% OUTPUT (struct). For each requested norm X in {L1, L2}:
%   bmr.X.kcal   BMR (kcal/day), NaN if no usable solution
%   bmr.X.atp    total ATP hydrolysed (mmol/day), NaN if none
%   bmr.X.stat   solver status (1 optimal; 3 suboptimal/numeric; -1/0 no solution)
%   bmr.X.v      flux vector (empty if none)
% plus bmr.atpMets : indices of the ATP metabolites summed.
%
% Anna Sheehy 2026 (companion to calculateBMR.m, which gives phenomenological
% reference BMRs from Harris-Benedict / Mifflin-St Jeor).

p = inputParser;
addParameter(p, 'brainATPfloor', 3500, @isnumeric);
addParameter(p, 'heartATPfloor', 1000, @isnumeric);
addParameter(p, 'fixObjective',  true, @islogical);
addParameter(p, 'kJperATP',      64,   @isnumeric);
addParameter(p, 'minNorm',       1e-6, @isnumeric);
addParameter(p, 'timeLimit',     1800, @isnumeric);
addParameter(p, 'lift',          true, @islogical);
addParameter(p, 'bigCoeff',      1e4,  @isnumeric);
addParameter(p, 'norms', {'L1','L2'}, @iscell);
parse(p, varargin{:});
opt = p.Results;

wantL1 = any(strcmpi(opt.norms, 'L1'));
wantL2 = any(strcmpi(opt.norms, 'L2'));

% --- build the BMR-determination model (as MethodSection3) -----------------
m = model;
m.lb(contains(m.rxns, 'Brain_DM_atp')) = opt.brainATPfloor;
m.lb(contains(m.rxns, 'Heart_DM_atp')) = opt.heartATPfloor;
if opt.fixObjective
    ob = strcmp(m.rxns, 'Whole_body_objective_rxn');
    m.lb(ob) = 1;  m.ub(ob) = 1;
end
m.c(:)      = 0;
m.osenseStr = 'min';

atpMets = find(contains(m.mets, '_atp['));     % all compartments; '_atp[' excludes datp
bmr     = struct('atpMets', atpMets);

% --- L1 (1-norm, LP) via optimizeCbModel -----------------------------------
% Sparse cross-check (under-reads); optimizeCbModel still honours the coupling.
if wantL1
    solL1  = optimizeCbModel(m, 'min', 'one');
    bmr.L1 = summariseBMR(solL1, m, atpMets, opt.kJperATP);
end

% --- L2 (2-norm, QP) — lift coupling + barrier-only, MethodSection3-faithful -
% The plain min-norm QP is numerically unstable on the WBM (Gurobi stat=-1 /
% CPLEX stat=3): badly-scaled COUPLING constraints wreck the barrier.
% liftCouplingConstraints (Sun, Fleming, Thiele 2013) reformulates them into a
% better-scaled, higher-dimensional problem; solving THAT barrier-only with
% crossover OFF gives stat=1 at the paper's minNorm=1e-6 for both sexes.
% (stat=3 is usually crossover failing to recover a vertex from an already-
% optimal barrier point; skipping it keeps the valid barrier solution.)
if wantL2
    m2 = m;
    if opt.lift && ~isModelLifted(m2) && isfield(m2, 'C') && ~isempty(m2.C)
        m2 = liftCouplingConstraints(m2, opt.bigCoeff, 0);   % printLevel 0; skip if already lifted
    end
    param2 = struct('minNorm', opt.minNorm, 'secondsTimeLimit', opt.timeLimit);
    if strcmpi(getCobraSolver('QP'), 'gurobi')
        param2.Method         = 2;     % barrier
        param2.Crossover      = 0;     % keep the barrier solution (no vertex recovery)
        param2.NumericFocus   = 3;     % maximum numerical care
        param2.BarHomogeneous = 1;     % homogeneous barrier: robust on hard problems
        param2.TimeLimit      = opt.timeLimit;
    end
    % Escalate minNorm only if a solve does not converge (SAFE: with c=0 the
    % min-norm solution is invariant to the weight; it only changes conditioning).
    solL2 = struct('stat', -1, 'v', []);
    for lam = opt.minNorm * [1 1e2 1e4]
        param2.minNorm = lam;
        solL2 = optimizeWBModel(m2, param2);
        if isfield(solL2, 'v') && ~isempty(solL2.v) && ismember(solL2.stat, [1 3])
            break
        end
    end
    bmr.L2 = summariseBMR(solL2, m, atpMets, opt.kJperATP);   % summariseBMR slices to original reactions
end
end


function r = summariseBMR(sol, m, atpMets, kJperATP)
% Convert a solved model into BMR by the Thiele 2020 §3.3 basal-metabolic-flux
% formula: consumption = [F R]*[vf; vr] with F=max(-S,0), R=max(S,0),
% vf=max(v,0), vr=max(-v,0); sum consumption over the ATP metabolites.
r = struct('stat', sol.stat, 'v', [], 'atp', NaN, 'kcal', NaN);
if ismember(sol.stat, [1 3]) && isfield(sol, 'v') && ~isempty(sol.v)
    v = sol.v(1:size(m.S,2));    % original reactions (lifted solves append aux variables)
    S = m.S;
    consumption = [max(-S, 0), max(S, 0)] * [max(v, 0); max(-v, 0)];   % §3.3
    r.v    = v;
    r.atp  = full(sum(consumption(atpMets)));                          % mmol ATP/day
    r.kcal = r.atp / 1000 * kJperATP * 0.239006;                      % 1 kJ = 0.239006 kcal
end
end
