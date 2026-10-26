function map = metabBoundMap(model, ipar, mets)
% metabBoundMap  Every bound that metabolomic personalisation of blood sets, as value = k * c.
%
%   map = metabBoundMap(model, ipar, mets)
%
% Calls the unchanged engine used by persWBMmetabolomics (Persephone),
% physiologicalConstraintsHMDBbased(..., 'direct', ..., 'bc'), twice with measured
% (min, max) = (1, 2) and (1, 3) uM for every metabolite. A bound that differs between
% the two calls depends on the maximum concentration (k = valueA/2); otherwise on the
% minimum (k = valueA). Any individual's bounds then follow from applyMetabBounds.
%
% INPUT   model  WBM;  ipar  IndividualParameters used for the model;  mets  VMH ids
% OUTPUT  map    table: rxnIdx, rxn, bound ('lb'/'ub'), dependsOn ('min'/'max'), k,
%                defaultValue, met, type (organ uptake / kidney filtration / blood-brain barrier)
%
% Anna Sheehy 2026.

n = numel(mets);
dA = [{'met','minCon','maxCon'}; [mets(:), num2cell(ones(n,1)), num2cell(2*ones(n,1))]];
dB = [{'met','minCon','maxCon'}; [mets(:), num2cell(ones(n,1)), num2cell(3*ones(n,1))]];
mA = physiologicalConstraintsHMDBbased(model, ipar, '', 'direct', dA, 'bc');
mB = physiologicalConstraintsHMDBbased(model, ipar, '', 'direct', dB, 'bc');
rows = cell(0, 7);
for side = {'lb', 'ub'}
    a = mA.(side{1});  b = mB.(side{1});  o = model.(side{1});
    for i = find(a ~= o | b ~= o)'
        if b(i) ~= a(i), k = a(i) / 2; dep = 'max'; else, k = a(i); dep = 'min'; end
        rows(end+1, :) = {i, model.rxns{i}, side{1}, dep, k, o(i), attributeMet(model, i, mets)}; %#ok<AGROW>
    end
end
map = cell2table(rows, 'VariableNames', {'rxnIdx','rxn','bound','dependsOn','k','defaultValue','met'});
noMet = cellfun(@isempty, map.met);
if any(noMet)
    warning('metabBoundMap:unattributed', '%d changed bounds not attributable to a panel metabolite, excluded: %s', ...
        sum(noMet), strjoin(unique(map.rxn(noMet)), ', '));
    map(noMet, :) = [];
end
map.bound = string(map.bound);
map.type = repmat({'organ uptake'}, height(map), 1);
map.type(startsWith(map.rxn, 'Kidney_')) = {'kidney filtration'};
map.type(startsWith(map.rxn, 'BBB_')) = {'blood-brain barrier'};
end

function met = attributeMet(model, i, mets)
parts = model.mets(model.S(:, i) ~= 0);
met = '';
for q = 1:numel(mets)
    if any(strcmp(parts, [mets{q} '[bc]'])) || any(strcmp(parts, [mets{q} '[csf]']))
        met = mets{q};  return
    end
end
end
