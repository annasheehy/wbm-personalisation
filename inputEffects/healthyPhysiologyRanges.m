function [lims, T] = healthyPhysiologyRanges(sourceFile, ipar)
% Low and high ends of the healthy range for each physiological input, in the units of
% IndividualParameters, from the values transcribed in sourceFile
% (references/healthyPhysiologyRanges.csv). Each row gives one field and sex as reported:
%   range    v1/v2 = low/high limit (2.5th/97.5th percentile or a defined healthy range)
%   meanSD   v1 = mean, v2 = SD, n; several rows of one field and sex (e.g. age bands)
%            are pooled, and the range is mean +/- 1.96 SD
%   meanP5P95  v1 = mean, v2/v3 = 5th/95th percentile, n; SD = (v3 - v2) / (2 * 1.645),
%            then pooled as meanSD (normal distribution assumed)
%   lognP10P90  v1/v3 = 10th/90th percentile (v2 = median, recorded only), n; a log-normal
%            per row through its 10th and 90th percentile, rows pooled by n, and the range
%            is the 10th to 90th percentile of the pooled distribution (the reference
%            range of the source)
% Units converted to the model's: kg/m2 (BMI) -> kg at the reference height;
% mL/min/1.73m2 -> mL/min with the Du Bois body surface area of the reference individual.
%
% INPUT:  ipar = struct('male', IndividualParameters, 'female', IndividualParameters)
% OUTPUT: lims.(sex).(field) = [low high]; field HeightBSA repeats Height.
%         T = one row per field and sex with the limits before and after conversion.
%
% Anna Sheehy 2026.

z = 1.959964;
S = readtable(sourceFile, 'TextType', 'char', 'Delimiter', ',');
[g, field, sex] = findgroups(S.field, S.sex);
lims = struct('male', struct(), 'female', struct());
T = table('Size', [0 7], 'VariableTypes', {'string','string','double','double','string','double','double'}, ...
    'VariableNames', {'field','sex','lowSource','highSource','unit','low','high'});
for k = 1:max(g)
    R = S(g == k, :);
    switch R.stat{1}
        case 'range'
            lim = [R.v1(1), R.v2(1)];
        case {'meanSD', 'meanP5P95'}                    % pool rows (e.g. age bands) by n
            sdRow = R.v2;
            if strcmp(R.stat{1}, 'meanP5P95'), sdRow = (R.v3 - R.v2) / (2 * 1.644854); end
            N = sum(R.n);  m = sum(R.n .* R.v1) / N;
            sd = sqrt(sum(R.n .* (sdRow.^2 + (R.v1 - m).^2)) / N);
            lim = m + [-z, z] * sd;
        case 'lognP10P90'
            mu = (log(R.v1) + log(R.v3)) / 2;  sg = (log(R.v3) - log(R.v1)) / (2 * 1.281552);  w = R.n / sum(R.n);
            cdf = @(x) sum(w .* 0.5 .* erfc(-(log(x) - mu) ./ (sg * sqrt(2))));
            lim = [fzero(@(x) cdf(x) - 0.10, median(R.v1)), fzero(@(x) cdf(x) - 0.90, median(R.v3))];
        otherwise
            error('healthyPhysiologyRanges: unknown stat ''%s'' for %s', R.stat{1}, field{k});
    end
    ip = ipar.(sex{k});
    switch R.unit{1}
        case 'kg/m2',          val = lim * (ip.Height / 100)^2;
        case 'mL/min/1.73m2',  val = lim * bsaDuBois(ip.bodyWeight, ip.Height) / 1.73;
        otherwise,             val = lim;
    end
    lims.(sex{k}).(field{k}) = val;
    T(end+1, :) = {string(field{k}), string(sex{k}), lim(1), lim(2), string(R.unit{1}), val(1), val(2)}; %#ok<AGROW>
end
for s = {'male', 'female'}
    if isfield(lims.(s{1}), 'Height'), lims.(s{1}).HeightBSA = lims.(s{1}).Height; end
end
end
