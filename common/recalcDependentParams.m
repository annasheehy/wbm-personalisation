function ip = recalcDependentParams(ip, field, val)
% Update one primary field in an IndividualParameters struct and
% recompute all downstream dependent parameters.
%
% INPUT
%   ip     IndividualParameters struct (from standardPhysiolDefaultParameters)
%   field  Name of the field to change (string)
%   val    New value for that field
%
% OUTPUT
%   ip     Updated struct with field set to val and dependents recalculated

ip.(field)       = val;
ip.CardiacOutput = ip.HeartRate * ip.StrokeVolume;

if strcmp(ip.sex, 'male')
    ip.BloodVolume = (0.3669*(ip.Height/100)^3 + 0.03219*ip.bodyWeight + 0.6041)*1000;
else
    ip.BloodVolume = (0.3561*(ip.Height/100)^3 + 0.03308*ip.bodyWeight + 0.1833)*1000;
end
end
