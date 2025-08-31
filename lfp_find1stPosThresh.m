function result = lfp_find1stPosThresh(filenum, thresh)
%lfp_find1stThresh for lfp_createEvents finds the first positive threshold
%   crossing after the alignment event in each trial.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

lfp_declareGlobals;

xingidx = zeros(size(lfp_SelectedTrials));
for trial = 1:length(lfp_SelectedTrials)
    refidx = find(ismember( ...
        lfp_Events(lfp_TrialIndex(trial,1):lfp_TrialIndex(trial,2), 2), ...
        lfp_AlignmentRef )) + lfp_TrialIndex(trial,1) - 1;
    if isempty(refidx)
        warning('lfp_find1stPosThresh:noref', ...
            'Trial %d has no reference event.', trial);
    else
        reftime = lfp_Events(refidx(1),1);
        refsample = lfp_time2index(reftime);
        xings = find( ...
            (lfp_Samples{filenum}(refsample + 1 : lfp_TrialIndex(trial,4)) > thresh) ...
            & (lfp_Samples{filenum}(refsample : lfp_TrialIndex(trial,4) - 1) < thresh) ...
            ) + refsample - 1;
        xingidx(trial) = xings(1);
    end
end
xingidx(xingidx==0) = [];
result = lfp_index2time(xingidx) + lfp_SamplePeriod/2;
