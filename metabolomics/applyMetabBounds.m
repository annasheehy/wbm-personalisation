function m = applyMetabBounds(model, map, metIDs, conc, w)
% applyMetabBounds  Bounds of one metabolomics-personalised WBM from a metabBoundMap.
%
%   m = applyMetabBounds(model, map, metIDs, conc, w)
%
% conc (uM) aligned with metIDs; w = measurement window (persWBMmetabolomics uses 0.2,
% i.e. [0.8c, 1.2c]). Bound = k * c * (1+w) where it depends on the maximum
% concentration, k * c * (1-w) where it depends on the minimum. Rows of map whose
% metabolite is not in metIDs are ignored (that metabolite stays at its default).
%
% Anna Sheehy 2026.

[inSet, loc] = ismember(map.met, metIDs);
map = map(inSet, :);  loc = loc(inSet);
c = conc(loc);  c = c(:);
f = (1 + w) * strcmp(map.dependsOn, 'max') + (1 - w) * strcmp(map.dependsOn, 'min');
val = map.k .* c .* f;
m = model;
isLb = map.bound == "lb";
m.lb(map.rxnIdx(isLb))  = val(isLb);
m.ub(map.rxnIdx(~isLb)) = val(~isLb);
end
