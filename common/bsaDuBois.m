function bsa = bsaDuBois(bodyWeight, height)
% Body surface area (m2) by the Du Bois and Du Bois (1916) formula, as used by
% individualizedLabReport (cardiac output option 4).
%
% INPUT:  bodyWeight  kg;  height  cm
% OUTPUT: bsa         m2
%
% Anna Sheehy 2026.

bsa = 0.007184 * bodyWeight.^0.425 .* height.^0.725;
end
