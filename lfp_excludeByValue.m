function lfp_excludeByValue(trials, filenums, upperlim, lowerlim)
%lfp_excludeByValue(trials, filenums, upperlim, lowerlim)
% Adds to lfp_BadTrials any trials that contain data points with
% values >= upperlim or <= lowerlim within lfp_findCommonTime(trials).

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
[interval, rawtrialinfo] = lfp_findCommonTime(trials);
xlimpoints = round(lfp_XLimAll/lfp_SamplePeriod);
interval(1) = max(interval(1), xlimpoints(1));
interval(2) = min(interval(2), xlimpoints(2));
if interval(1) >= interval(2)
    error('lfp_excludeByValue:nodata', ...
        'There is no common time among trials and lfp_XLimAll.' );
end
for trialidx = 1:length(trials)
    trial = trials(trialidx);
    refpoint = rawtrialinfo(trialidx,3);
    if refpoint == 0
        warning('lfp_excludeByValue:skipping', ...
            'Skipping trial %d, no ref event', trial );
        continue
    end
    samplerange = refpoint + (interval(1):interval(2));
    for filenum = filenums
        outofbounds = find( ...
            lfp_Samples{filenum}(samplerange) >= upperlim ...
            | lfp_Samples{filenum}(samplerange) <= lowerlim );
        if ~isempty(outofbounds)
            lfp_BadTrials(end+1) = trial;
            break
        end
    end
end

