function model = buildFedWBM(rawModel, IndividualParameters, Diet, organWeights, hctRef)
% Fed whole-body model from a raw Harvey/Harvetta and a set of IndividualParameters, built
% as the control models of fastingWBM.m (Parts 1-2): physiological constraints
% (physiologicalConstraintsHMDBbased), whole-body maintenance as objective, organ stores
% set to fill (setFeedingFastingConstraints, feeding) and the diet applied
% (setDietConstraints). With the reference parameters this reproduces the control bounds.
%
% INPUT:  rawModel              Harvey/Harvetta as loaded (liver biomass closed, fastingWBM.m)
%         IndividualParameters  as from standardPhysiolDefaultParameters, personalised
%         Diet                  diet cell array (e.g. EUAverageDietNew)
%         organWeights          optional: 'none' (default) keeps the organ weights of the
%                               model; 'std' recalculates them from the parameters
%                               (calcOrganFract, as persWBM); 'pers' does so with the
%                               personalised blood volume (calcOrganFractPersonal)
%         hctRef                optional, for 'pers': reference haematocrit to which the red
%                               blood cell weight is scaled
% OUTPUT: model                 fed model (not lifted; see onLiftedModel)
%
% Anna Sheehy 2026.

if nargin < 4 || isempty(organWeights), organWeights = 'none'; end
if nargin < 5, hctRef = []; end
model = rawModel;
switch organWeights
    case 'none'
    case 'std'
        [organs, ~, OrganWeightFract, IndividualParameters] = calcOrganFract(model, IndividualParameters);
        model = adjustWholeBodyRxnCoeff(model, organs, OrganWeightFract);
    case 'pers'
        [organs, ~, OrganWeightFract, IndividualParameters] = calcOrganFractPersonal(model, IndividualParameters, hctRef);
        model = adjustWholeBodyRxnCoeff(model, organs, OrganWeightFract);
    otherwise
        error('buildFedWBM: organWeights must be ''none'', ''std'' or ''pers''');
end
model = physiologicalConstraintsHMDBbased(model, IndividualParameters);
model = changeObjective(model, 'Whole_body_objective_rxn');
model = setDietConstraints(setFeedingFastingConstraints(model, 'feeding', -1e6, 1e6), Diet);
model.osenseStr = 'max';
end
