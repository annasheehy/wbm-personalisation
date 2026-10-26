function [organs, OrganWeight, OrganWeightFract, ip] = calcOrganFractPersonal(model, ip, hctRef)
% Organ weights for the whole-body maintenance reaction, as calcOrganFract (PSCM toolbox), but
% with the blood cell weights calculated from the personalised blood volume. calcOrganFract
% replaces ip.BloodVolume (e.g. from height and weight by the Nadler equation, or measured)
% with an estimate from body weight alone; this function keeps the personalised value.
% calcOrganFract itself is not changed, so persWBM (Persephone) behaves as published.
%
% INPUT:  model   whole-body model; ip  IndividualParameters with BloodVolume (mL)
%         hctRef  optional: if given, the red blood cell weight is also scaled by
%                 ip.Hematocrit / hctRef (red cell mass follows haematocrit); omit to keep
%                 the calcOrganFract assumption of 495 g red cells per L of blood.
% OUTPUT: as calcOrganFract; ip.BloodVolume is the personalised value used.
%
% Blood cell weights follow calcOrganFract: white cells 1% of blood weight (lymphocytes 30%
% of white cells; B cells 9% and NK cells 15% of lymphocytes; CD4 T cells 15% and monocytes
% 5% of white cells), platelets 4 g/L and red cells 495 g/L of blood, blood density 1.0506.
%
% Anna Sheehy 2026.

bv = ip.BloodVolume;                                      % mL, personalised
[organs, OrganWeight, OrganWeightFract, ip] = calcOrganFract(model, ip);
Wt = ip.bodyWeight * 1000;                                % g
bloodWeight = bv * 1.0506;                                % g
wbc = 0.01 * bloodWeight;  lymph = 0.3 * wbc;
rbcPerL = 495;
if nargin > 2 && ~isempty(hctRef), rbcPerL = 495 * ip.Hematocrit / hctRef; end
w = struct('Blood', bloodWeight, 'Bcells', 0.09 * lymph, 'CD4Tcells', 0.15 * wbc, ...
    'Nkcells', 0.15 * lymph, 'Monocyte', 0.05 * wbc, 'Platelet', 4 * bv / 1000, 'RBC', rbcPerL * bv / 1000);
for f = fieldnames(w)'
    i = find(strcmp(organs, f{1}));
    if ~isempty(i)
        OrganWeight(i, 1) = w.(f{1});
        OrganWeightFract(i, 1) = w.(f{1}) / Wt;
        ip.OrgansWeights(i, :) = {organs{i}, num2str(OrganWeight(i)), num2str(OrganWeightFract(i) * 100)};
    end
end
ip.BloodVolume = bv;
end
