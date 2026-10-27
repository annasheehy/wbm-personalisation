function caps = euAverageIonCaps()
% Uptake limits {rxn, lb} for water and the 11 non-phosphate dietary ions at their provision
% in the average European diet (mmol/day, from EUAverageDietNew). An alternative to the
% default of setFastingConstraints, which leaves them freely available:
%   setFastingConstraints(model, sex, 'ionCaps', euAverageIonCaps())
%
% Anna Sheehy 2026.

caps = {'Diet_EX_h2o[d]', -165892.3425;  'Diet_EX_na1[d]', -172.4463;  'Diet_EX_k[d]',  -83.4338; ...
        'Diet_EX_cl[d]',    -187.7787;   'Diet_EX_so4[d]',        0;    'Diet_EX_zn2[d]',  -0.1510; ...
        'Diet_EX_ca2[d]',    -18.5491;   'Diet_EX_cu2[d]',   -0.0357;   'Diet_EX_mn2[d]',  -0.0747; ...
        'Diet_EX_fe2[d]',     -0.1092;   'Diet_EX_fe3[d]',   -0.1092};
end
