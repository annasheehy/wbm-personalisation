function [kcal, atp] = bmrFromFlux(model, v, kJperATP)
% BMR of a flux distribution by the basal metabolic flux formula of Thiele 2020 (§3.3):
% the ATP consumed by every reaction in every organ and compartment is summed and
% converted to kcal. Consumption = [F R] * [vf; vr] with F = max(-S, 0), R = max(S, 0),
% vf = max(v, 0), vr = max(-v, 0).
%
% INPUT:  model     whole-body model;  v  flux vector (extra entries of a lifted model are ignored)
%         kJperATP  optional, energy per mol ATP hydrolysed (default 64, Wackerhage 1998)
% OUTPUT: kcal      kcal/day;  atp  ATP consumed (mmol/day)
%
% Anna Sheehy 2026.

if nargin < 3, kJperATP = 64; end
v = v(1:size(model.S, 2));
v = v(:);
atpMets = contains(model.mets, '_atp[');                  % all organs and compartments; excludes datp
consumption = [max(-model.S(atpMets, :), 0), max(model.S(atpMets, :), 0)] * [max(v, 0); max(-v, 0)];
atp  = full(sum(consumption));                            % mmol ATP/day
kcal = atp / 1000 * kJperATP * 0.239006;                  % mol x kJ/mol x kcal/kJ
end
