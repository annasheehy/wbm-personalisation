function setupPersWBM(projectFolder)
% Adds the persWBM folders to the MATLAB path. Run once per session, after
% initCobraToolbox. The archive folder is not added.
%
% INPUT: projectFolder (optional)  folder for the input tables, results and manuscript
%        files (see persWBMpaths); stored in the environment variable PERSWBM_PROJECT
%
% Anna Sheehy 2026.

here = fileparts(mfilename('fullpath'));
for d = {'common', 'inputEffects', 'metabolomics', 'fasting', 'active'}
    addpath(fullfile(here, d{1}));
end
if nargin > 0, setenv('PERSWBM_PROJECT', projectFolder); end
end
