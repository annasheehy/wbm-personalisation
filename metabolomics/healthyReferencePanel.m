function panel = healthyReferencePanel(sourceFile)
% Healthy-adult blood concentrations as log-normal distributions per sex, from the values
% transcribed in sourceFile (references/healthyReferenceValues.csv). Each row gives one
% metabolite for 'male', 'female' or 'both', in the form its source reports it:
%   meanSD   v1 = mean, v2 = SD                          moment-matched log-normal
%   meanCI   v1 = mean, v2/v3 = 95% CI of the mean, n    SD = sqrt(n) * (v3 - v2) / 3.92
%   pctl     v1 = 2.5th, v2 = 97.5th percentile          log-normal through both
%
% Output: one row per metabolite: vmh, name, hmdb, source, sexSpecific, and the log-scale
% parameters (concentrations in uM) muMale, sigmaMale, muFemale, sigmaFemale, with the
% 2.5th and 97.5th percentiles they imply (p025Male, p975Male, p025Female, p975Female).
%
% Anna Sheehy 2026.

z = 1.959964;
T = readtable(sourceFile, 'TextType', 'char', 'Delimiter', ',');
toUM = ones(height(T), 1);
toUM(strcmp(T.unit, 'mmol/L')) = 1000;
if ~all(ismember(T.unit, {'umol/L', 'mmol/L'}))
    error('healthyReferencePanel: unknown unit in %s', sourceFile);
end
v1 = T.v1 .* toUM;  v2 = T.v2 .* toUM;  v3 = T.v3 .* toUM;

mu = nan(height(T), 1);  sigma = nan(height(T), 1);
for r = 1:height(T)
    switch T.stat{r}
        case 'meanSD'
            m = v1(r);  sd = v2(r);
        case 'meanCI'
            m = v1(r);  sd = sqrt(T.n(r)) * (v3(r) - v2(r)) / (2 * z);
        case 'pctl'
            mu(r) = (log(v1(r)) + log(v2(r))) / 2;
            sigma(r) = (log(v2(r)) - log(v1(r))) / (2 * z);
            continue
        otherwise
            error('healthyReferencePanel: unknown stat ''%s'' for %s', T.stat{r}, T.vmh{r});
    end
    sigma(r) = sqrt(log(1 + (sd / m)^2));
    mu(r) = log(m) - sigma(r)^2 / 2;
end

[vmh, first] = unique(T.vmh, 'stable');
n = numel(vmh);
panel = table(vmh, T.name(first), T.hmdb(first), T.source(first), false(n, 1), ...
    nan(n, 1), nan(n, 1), nan(n, 1), nan(n, 1), 'VariableNames', ...
    {'vmh', 'name', 'hmdb', 'source', 'sexSpecific', 'muMale', 'sigmaMale', 'muFemale', 'sigmaFemale'});
for q = 1:n
    rows = find(strcmp(T.vmh, vmh{q}));
    for r = rows(:)'
        switch T.sex{r}
            case 'both',   panel.muMale(q) = mu(r);  panel.sigmaMale(q) = sigma(r);
                           panel.muFemale(q) = mu(r);  panel.sigmaFemale(q) = sigma(r);
            case 'male',   panel.muMale(q) = mu(r);  panel.sigmaMale(q) = sigma(r);
            case 'female', panel.muFemale(q) = mu(r);  panel.sigmaFemale(q) = sigma(r);
            otherwise, error('healthyReferencePanel: unknown sex ''%s'' for %s', T.sex{r}, vmh{q});
        end
    end
    panel.sexSpecific(q) = ~any(strcmp(T.sex(rows), 'both'));
end
if any(isnan([panel.muMale; panel.muFemale]))
    error('healthyReferencePanel: missing male or female values for %s', ...
        strjoin(panel.vmh(isnan(panel.muMale) | isnan(panel.muFemale)), ', '));
end
panel.p025Male   = exp(panel.muMale - z * panel.sigmaMale);
panel.p975Male   = exp(panel.muMale + z * panel.sigmaMale);
panel.p025Female = exp(panel.muFemale - z * panel.sigmaFemale);
panel.p975Female = exp(panel.muFemale + z * panel.sigmaFemale);
end
