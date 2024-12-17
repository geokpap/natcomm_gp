function [startrecsample, endrecsample, recOnSamples, recOffSamples] ...
    = lfp_findRecSegments
% Find the beginning and end of the recorded time interval that contains
% each trial, using the values in lfp_TimeStamps.  For this purpose,
% "contains the trial" is defined as "contains the trial's
% lfp_NominalTrialStart event".  <startrecsample> and <endrecsample> each
% have one element for each trial.  <recOnSamples> and <recOffSamples> each
% have one element for each gap found plus one for start or end of session.

%$Rev: 369 $
%$Date: 2015-10-15 19:09:39 -0400 (Thu, 15 Oct 2015) $
%$Author: dgibson $

global lfp_NoWaitbar lfp_SamplesPerFrame lfp_TimeStamps lfp_TrialIndex

if ~lfp_NoWaitbar
    hWaitBar = waitbar(0, '', 'Name', 'Finding recorded segments');
    WaitBarSteps = size(lfp_TrialIndex,1);
end

% Start by finding any gaps in the recorded sample data, where a gap is
% defined as any time interval between successive timestamps that differs
% from the normal frame rate by more than half a frame.  The "normal" frame
% rate is estimated as the median.
difs = reshape(lfp_TimeStamps(2:end) - lfp_TimeStamps(1:end-1), ...
    1, []);
meddif = median(difs);
gaps = find((difs < meddif/2) | (difs > meddif*1.5));
if ~lfp_NoWaitbar
    waitbar(1/WaitBarSteps, hWaitBar);
end

% Then convert to the sample indices of the first and last samples in each
% recorded segment, based on the assumption that an entire frame is
% recorded at a time (which seems to be true looking at the data in Dec
% 2003). 
numsamples = length(lfp_TimeStamps) * lfp_SamplesPerFrame;
recOnSamples = [0 gaps] * lfp_SamplesPerFrame + 1;
recOffSamples = [gaps * lfp_SamplesPerFrame, numsamples];

% Now we can find the last Record On before the start of each trial and the
% first Record Off after the end of each trial. 
startrecsample = zeros(size(lfp_TrialIndex,1), 1);
endrecsample = zeros(size(lfp_TrialIndex,1), 1);
for trialnum = 1:size(lfp_TrialIndex,1)
    if ~lfp_NoWaitbar
        waitbar(trialnum/WaitBarSteps, hWaitBar);
    end
    found = recOnSamples(recOnSamples ...
        <= lfp_TrialIndex(trialnum, 3) );
    startrecsample(trialnum) = found(end);
    found = recOffSamples(recOffSamples ...
        >= lfp_TrialIndex(trialnum, 4) );
    endrecsample(trialnum) = found(1);
end
if ~lfp_NoWaitbar
    close(hWaitBar);
end
recOnSamples = reshape(recOnSamples, [], 1);
recOffSamples = reshape(recOffSamples, [], 1);
