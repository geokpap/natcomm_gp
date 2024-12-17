function [startSampleIndex, endSampleIndex] = ...
    lfp_addTIcols3n4(trialstartindex, trialendindex, events);
% <trialstartindex> and <trialendindex> are vectors of lfp_Event row
% numbers, no special orientation required.
% Find the corresponding sample indices.  The convention is
% to use the CSC timestamps at or immediately before the start event, and
% at or immediately after the end event (if this would  run off the end of
% the data, then the last sample is used).  For further speed
% improvement, I rely on the members of trialstartindex and trialendindex
% having nondecreasing timestamps when you look at successive start, end
% pairs.
% [startSampleIndex, endSampleIndex] are column vectors.

%$Rev: 56 $
%$Date: 2009-04-10 13:51:07 -0400 (Fri, 10 Apr 2009) $
%$Author: dgibson $

lfp_declareGlobals;

if ~lfp_NoWaitbar
    hWaitBar = waitbar(0, '', 'Name', 'Constructing Trial Index');
    WaitBarSteps = length(trialstartindex);
    CurrentWaitBarStep = 0;
end
startSampleIndex = zeros(length(trialstartindex), 1);
endSampleIndex = zeros(length(trialstartindex), 1);
previousindex = 1;
previousTSindex = 1;
totalsamples = lfp_SamplesPerFrame * length(lfp_TimeStamps);
for k = 1:length(trialstartindex)
    startSampleIndex(k) = lfp_time2index(events(trialstartindex(k),1));
    previousindex = trialstartindex(k);
    previousTSindex = startSampleIndex(k);
    endSampleIndex(k) = lfp_time2index(events(trialendindex(k),1));
    if endSampleIndex(k) > totalsamples
        endSampleIndex(k) = totalsamples;
    elseif lfp_index2time(endSampleIndex(k) - 1) == events(trialendindex(k),1)
        endSampleIndex(k) = endSampleIndex(k) - 1;
    end
    previousindex = trialendindex(k);
    previousTSindex = endSampleIndex(k);
    if ~lfp_NoWaitbar
        CurrentWaitBarStep = CurrentWaitBarStep + 1;
        waitbar(CurrentWaitBarStep/WaitBarSteps, hWaitBar);
    end
end
if ~lfp_NoWaitbar
    close(hWaitBar);
end
