function model = onLiftedModel(liftedModel, newModel)
% Puts the bounds, objective and changed stoichiometric columns of newModel onto
% liftedModel, a model with the same reactions whose coupling constraints have already been
% reformulated by liftCouplingConstraints. Lifting changes only the coupling constraints,
% so the result equals lifting newModel, without the slow lift.
%
% INPUT:  liftedModel  lifted model (e.g. a control model of fastingWBM.m)
%         newModel     the same model rebuilt with other parameters (e.g. buildFedWBM)
% OUTPUT: model        lifted model with the bounds of newModel
%
% Only the whole-body maintenance column changes when organ weights change, so at most
% one column may differ; anything else is an error.
%
% Anna Sheehy 2026.

assert(isequal(liftedModel.rxns, newModel.rxns) && isequal(liftedModel.mets, newModel.mets), ...
    'onLiftedModel: reaction or metabolite lists differ');
model = liftedModel;
model.lb = newModel.lb;  model.ub = newModel.ub;  model.c = newModel.c;
d = find(any(model.S ~= newModel.S, 1));
if numel(d) > 1
    error('onLiftedModel: %d stoichiometric columns differ (expected 0 or 1): %s', ...
        numel(d), strjoin(model.rxns(d), ', '));
end
model.S(:, d) = newModel.S(:, d);
end
