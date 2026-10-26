function tf = isModelLifted(model)
% isModelLifted  True if a model has already been through liftCouplingConstraints.
%
%   tf = isModelLifted(model)
%
% Use before calling liftCouplingConstraints so a model is never lifted twice
% (lifting is the expensive step, and it caches into the pipeline checkpoints).
% Markers: the added variables/matrix (evars/E) or a '_lifted' modelID stamp.
%
% Anna Sheehy 2026.
tf = (isfield(model, 'evars') && ~isempty(model.evars)) || ...
     (isfield(model, 'E')     && ~isempty(model.E))     || ...
     (isfield(model, 'modelID') && contains(lower(char(string(model.modelID))), 'lifted'));
end
