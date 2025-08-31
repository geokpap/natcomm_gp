function [samplerange, refpoint] = lfp_getSampleRange(trial)
%samplerange = lfp_getSampleRange(trial)
%   Returns the range of absolute sample indices that are within the
%   recorded segment containing <trial>, and that are also within
%   lfp_XLimAll if lfp_XLimAll is not empty, or within
%   lfp_NominalTrialStart and lfp_NominalTrialEnd if lfp_XLimAll is empty.
%   Returns [] if there is no lfp_AlignmentRef in <trial>.  Also returns
%   the absolute index of the sample closest to the lfp_AlignmentRef.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

lfp_declareGlobals;

trialevents = lfp_Events( ...
    lfp_TrialIndex(trial,1) : lfp_TrialIndex(trial,2), ...
    : );
reftime = trialevents( ...
    find(ismember(trialevents(:,2), lfp_AlignmentRef)), ...
    1 );
if length(reftime) == 0
    samplerange = [];
    refpoint = 0;
    return
else
    refpoint = lfp_time2index(reftime(1));
end
if isempty(lfp_XLimAll)
    endpoints = lfp_TrialIndex(trial,3:4);
else
    xlimpoints = round(lfp_XLimAll/lfp_SamplePeriod);
    endpoints = xlimpoints + refpoint;
end
endpoints(1) = max(endpoints(1), lfp_TrialRec(trial,1));
endpoints(2) = min(endpoints(2), lfp_TrialRec(trial,2));
samplerange = endpoints(1) : endpoints(2);
    
