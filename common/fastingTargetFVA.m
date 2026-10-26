function F = fastingTargetFVA(model, rxns, varargin)
% fastingTargetFVA  Flux variability (min / max) of selected reactions under BMR conditions.
%
%   F = fastingTargetFVA(model, rxns)
%   F = fastingTargetFVA(model, rxns, 'Name', Value, ...)
%
% For each reaction in rxns, solves two LPs (minimise, then maximise its flux) on the
% SAME constrained model that calculateModelBMR uses for the L2 (min-norm) solution:
% resting brain and heart ATP floors, whole-body objective fixed to 1. The result is
% the range of flux the reaction CAN carry while every other constraint (network,
% coupling, diet/storage, floors) holds - i.e. whether a literature value is
% reachable at all - which the nominal lb/ub alone cannot tell you.
%
% The L2 point from calculateModelBMR is the flux the model actually CHOSES; this is
% the range it could have chosen from. Plot them together.
%
% INPUT
%   model   whole-body model (e.g. the fasted Harvey/Harvetta from setFastingConstraints)
%   rxns    cell/string array of reaction ids
%
% Name-Value options (defaults identical to calculateModelBMR):
%   'brainATPfloor'  3500    lb on reactions matching 'Brain_DM_atp'
%   'heartATPfloor'  1000    lb on reactions matching 'Heart_DM_atp'
%   'fixObjective'   true    Whole_body_objective_rxn fixed to 1
%
% OUTPUT F (table): rxn, fvaMin, fvaMax, statMin, statMax, seconds
%   fvaMin/fvaMax are NaN where a solve did not return an optimal solution (stat ~= 1),
%   or where the reaction is not in the model (warned, never silently skipped).
%
% Anna Sheehy 2026. Companion to calculateModelBMR.m / fastingFigures.m.

p = inputParser;
addParameter(p, 'brainATPfloor', 3500, @isnumeric);
addParameter(p, 'heartATPfloor', 1000, @isnumeric);
addParameter(p, 'fixObjective',  true, @islogical);
parse(p, varargin{:});
opt = p.Results;

% --- the BMR-determination model, exactly as calculateModelBMR builds it ----------
m = model;
m.lb(contains(m.rxns, 'Brain_DM_atp')) = opt.brainATPfloor;
m.lb(contains(m.rxns, 'Heart_DM_atp')) = opt.heartATPfloor;
if opt.fixObjective
    ob = strcmp(m.rxns, 'Whole_body_objective_rxn');
    m.lb(ob) = 1;  m.ub(ob) = 1;
end

rxns = cellstr(rxns);
n = numel(rxns);
F = table(string(rxns(:)), nan(n,1), nan(n,1), nan(n,1), nan(n,1), nan(n,1), ...
    'VariableNames', {'rxn','fvaMin','fvaMax','statMin','statMax','seconds'});

for k = 1:n
    idx = find(strcmp(m.rxns, rxns{k}), 1);
    if isempty(idx)
        warning('fastingTargetFVA:missing', 'Reaction "%s" not in model - skipped. Check against model.rxns.', rxns{k});
        continue
    end
    t0 = tic;
    m.c(:) = 0;  m.c(idx) = 1;
    solMin = optimizeCbModel(m, 'min');
    solMax = optimizeCbModel(m, 'max');
    F.statMin(k) = solMin.stat;  F.statMax(k) = solMax.stat;
    if solMin.stat == 1, F.fvaMin(k) = solMin.f; end
    if solMax.stat == 1, F.fvaMax(k) = solMax.f; end
    F.seconds(k) = toc(t0);
    fprintf('  FVA %-30s [%10.3f, %10.3f]  stat %d/%d  %.0f s\n', rxns{k}, F.fvaMin(k), F.fvaMax(k), ...
        solMin.stat, solMax.stat, F.seconds(k));
end
end
