function index = lfp_findLastPointByValue(...
    trials, filenums, upperlim, lowerlim, n )
%lfp_findLastPointByValue(trials, filenums, upperlim, lowerlim, n)
% Looks at each trial's wave data for each of <filenums> within
% lfp_findCommonTime(trials) to find the last data point with value >=
% <upperlim> or <= lowerlim, then searches trial by trial to find the <n>th
% latest such point, and returns the index to that point relative to
% lfp_AlignmentRef.  Returns -Inf if there fewer than <n> trials that have
% out-of-bounds data points.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

lfp_declareGlobals;

if isempty(trials)
    trials = lfp_enabledTrials(1:size(lfp_TrialIndex,1));
end
if isempty(filenums)
    filenums = lfp_ActiveFilenums;
end
lastpt(1:length(trials)) = -Inf;
[interval, rawtrialinfo] = lfp_findCommonTime(trials);
xlimpoints = round(lfp_XLimAll/lfp_SamplePeriod);
interval(1) = max(interval(1), xlimpoints(1));
interval(2) = min(interval(2), xlimpoints(2));
for trialidx = 1:length(trials)
    trial = trials(trialidx);
    refpoint = rawtrialinfo(trialidx,3);
    if refpoint == 0
        warning('lfp_findLastPointByValue:skipping', ...
            'Skipping trial %d, no ref event', trial );
        continue
    end
    samplerange = refpoint + (interval(1):interval(2));
    for filenum = filenums
        outofbounds = find( ...
            lfp_Samples{filenum}(samplerange) >= upperlim ...
            | lfp_Samples{filenum}(samplerange) <= lowerlim );
        if ~isempty(outofbounds)
            % make outofbounds relative to refpoint:
            outofbounds = outofbounds - (refpoint - samplerange(1) + 1);
            lastpt(trialidx) = max(lastpt(trialidx), outofbounds(end));
        end
    end
end
if length(lastpt) <  n
    index = -Inf;
else
    lastpt = sort(lastpt);
    index = lastpt(end - n + 1);
end

