%% oatBMRdecomposition.m
% Thesis Ch5: why heart rate and haematocrit move BMR in the one-at-a-time results.
% Decomposes the BMR change (QP-FBA, low -> high input value) by organ and by reaction,
% using the ATP consumption counted by calculateModelBMR, and prints whole-body O2 uptake,
% CO2 output and the brain and heart ATP demands. Reads the saved jobs of metabolomicsOAT.m
% (oatMode 'l2'); no solves.
%
% Anna Sheehy 2026.

P = persWBMpaths();
outDir = fullfile(P.results, 'metabolomics', 'oatL2');
inputs = [3 4];                                           % physInputs of metabolomicsOAT: heart rate, haematocrit
M = load(fullfile(P.results, 'fasting', 'ckpt_feasibleModels.mat'), 'maleFed', 'femaleFed');
mods = struct('male', M.maleFed, 'female', M.femaleFed);
kc = 64 * 0.239006 / 1000;                                % mmol ATP/day -> kcal/day
for s = {'male', 'female'}
    m = mods.(s{1});  n = numel(m.rxns);
    atp = contains(m.mets, '_atp[');
    Fa = full(sum(max(-m.S(atp, :), 0), 1))';  Ra = full(sum(max(m.S(atp, :), 0), 1))';
    cons = @(v) max(v, 0) .* Fa + max(-v, 0) .* Ra;          % ATP consumed by each reaction
    org = regexp(m.rxns, '^([A-Za-z0-9]+)_', 'tokens', 'once');
    org(cellfun(@isempty, org)) = {{'exchange'}};  org = cellfun(@(c) c{1}, org, 'UniformOutput', false);
    fprintf('\n===== %s\n', s{1});
    for id = inputs
        lo = load(fullfile(outDir, sprintf('%s_physiology_%02d_1.mat', s{1}, id)));
        hi = load(fullfile(outDir, sprintf('%s_physiology_%02d_2.mat', s{1}, id)));
        vL = double(lo.v(1:n));  vH = double(hi.v(1:n));
        dc = cons(vH) - cons(vL);
        fprintf('\n-- input %d: %.3g -> %.3g; BMR %.0f -> %.0f kcal/day\n', id, lo.val, hi.val, lo.kcal, hi.kcal);
        [u, ~, j] = unique(org);  byOrg = accumarray(j, dc) * kc;
        [~, o] = sort(abs(byOrg), 'descend');
        fprintf('   BMR change by organ (kcal/day): %s\n', strjoin(compose('%s %+.0f', string(u(o(1:8))), byOrg(o(1:8))), ', '));
        [~, o] = sort(abs(dc), 'descend');
        for r = o(1:12)'
            fprintf('   %-40s %10.1f -> %10.1f  ATP %+9.1f (%+.1f kcal)\n', m.rxns{r}, vL(r), vH(r), dc(r), kc * dc(r));
        end
        for x = {'EX_o2[a]', 'EX_co2[a]', 'Heart_DM_atp_c_', 'Brain_DM_atp_c_'}
            r = find(strcmp(m.rxns, x{1}));
            if ~isempty(r), fprintf('   %-20s %10.1f -> %10.1f\n', x{1}, vL(r), vH(r)); end
        end
    end
end
