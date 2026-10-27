function model = setBMRConditions(model, brainATPfloor, heartATPfloor)
% BMR conditions of the basal metabolic flux method (Thiele 2020; MethodSection3.mlx):
% whole-body maintenance fixed to one day (Whole_body_objective_rxn = 1) and minimum ATP
% hydrolysis in brain and heart (resting electrophysiology and cardiac work).
%
% INPUT:  model          whole-body model
%         brainATPfloor  optional, lb on Brain_DM_atp (mmol/day); default 3500, [] = unchanged
%         heartATPfloor  optional, lb on Heart_DM_atp (mmol/day); default 1000, [] = unchanged
% OUTPUT: model          with these bounds set
%
% Anna Sheehy 2026.

if nargin < 2, brainATPfloor = 3500; end
if nargin < 3, heartATPfloor = 1000; end
if ~isempty(brainATPfloor), model.lb(contains(model.rxns, 'Brain_DM_atp')) = brainATPfloor; end
if ~isempty(heartATPfloor), model.lb(contains(model.rxns, 'Heart_DM_atp')) = heartATPfloor; end
ob = strcmp(model.rxns, 'Whole_body_objective_rxn');
model.lb(ob) = 1;  model.ub(ob) = 1;
end
