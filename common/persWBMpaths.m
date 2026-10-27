function p = persWBMpaths()
% Folders used by the persWBM scripts. The project folder holds the literature input
% tables (references/), the results (results/) and the manuscript (ThesisChapter/LateX,
% main.tex and appendix_chapter5.tex), which receives the figures and the appendix. Generated
% LaTeX tables are written to results/latex, to be pasted into main.tex. It is taken from the environment variable
% PERSWBM_PROJECT (set it with setupPersWBM(folder)), or else from the author's machines.
%
% OUTPUT: p.project, p.references, p.results, p.latex, p.tables, p.figures
%
% Anna Sheehy 2026.

root = getenv('PERSWBM_PROJECT');
if isempty(root)
    known = {'D:\OneDrive - National University of Ireland, Galway\Group_MSP - Personalising WBMs', ...
             'C:\Users\AnnaS\OneDrive - National University of Ireland, Galway\Group_MSP - Personalising WBMs'};
    hit = find(cellfun(@(f) exist(f, 'dir') == 7, known), 1);
    if isempty(hit)
        error('persWBMpaths:noProject', 'Project folder not set: run setupPersWBM(folder) first.');
    end
    root = known{hit};
end
p.project    = root;
p.references = fullfile(root, 'references');
p.results    = fullfile(root, 'results');
p.latex      = fullfile(root, 'ThesisChapter', 'LateX');
p.tables     = fullfile(root, 'results', 'latex');
p.figures    = fullfile(p.latex, 'Figures');
end
