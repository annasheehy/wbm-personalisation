function R = oatReadouts()
% The 11 organ exchanges used as read-outs of the one-at-a-time analyses
% (metabolomicsOAT.m, heightCardiacOutput.m): the main fuels carried between organs.
% Liver and kidney exchanges are those with the circulating blood ([bc]); kidney glucose
% is glomerular filtration; brain exchanges are with the cerebrospinal fluid ([csf]).
%
% OUTPUT: R  11 x 2 cell, {reaction id, label}
%
% Anna Sheehy 2026.

R = {'Liver_EX_glc_D(e)_[bc]',     'Liver glucose'
     'Brain_EX_glc_D(e)_[csf]',    'Brain glucose'
     'Muscle_EX_glc_D(e)_[bc]',    'Muscle glucose'
     'RBC_EX_glc_D(e)_[bc]',       'RBC glucose'
     'Kidney_EX_glc_D(e)_[bc]',    'Kidney glucose'
     'Adipocytes_EX_glyc(e)_[bc]', 'Adipose glycerol'
     'Liver_EX_bhb(e)_[bc]',       'Liver 3-hydroxybutyrate'
     'Brain_EX_bhb(e)_[csf]',      'Brain 3-hydroxybutyrate'
     'Muscle_EX_ala_L(e)_[bc]',    'Muscle alanine'
     'Liver_EX_ala_L(e)_[bc]',     'Liver alanine'
     'EX_urea[u]',                 'Urinary urea'};
end
