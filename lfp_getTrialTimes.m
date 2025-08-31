function [reltrialstart, reltrialend, reftime] = lfp_getTrialTimes(trial, window)
%[reltrialstart, reltrialend, reftime] = lfp_getTrialTimes(trial)
%   Returns the start time and end time of trial number <trial> relative to
%   the lfp_AlignmentRef event, along with the absolute timestamp of the
%   lfp_AlignmentRef event.  If <window> is not empty, it contains start and
%   end times relative to the lfp_AlignmentRef event, and those values
%   override the actual start and end times.  Trial times are constrained
%   to stay within the corresponding lfp_TrialRec bounds.

%$Rev: 351 $
%$Date: 2015-05-28 12:32:51 -0400 (Thu, 28 May 2015) $
%$Author: dgibson $

global lfp_TrialIndex lfp_Events lfp_AlignmentRef lfp_TrialRec

eventrange = lfp_TrialIndex(trial,1) : lfp_TrialIndex(trial,2);
trialevents = lfp_Events(eventrange,:);
reftime = trialevents( ...
    find(ismember(trialevents(:,2), lfp_AlignmentRef)), ...
    1 );
if isempty(reftime)
    error('lfp_getTrialTimes:noRefEvent', ...
        [ 'Could not find the reference event in trial ' ...
            num2str(trial) ]);
else
    reftime = reftime(1);
end
if isempty(window)
    reltrialstart = lfp_Events(lfp_TrialIndex(trial,1), 1) - reftime;
    reltrialend = lfp_Events(lfp_TrialIndex(trial,2), 1) - reftime;
else
    reltrialstart = window(1);
    reltrialend = window(2);
end
reltrialstart = max(reltrialstart, ...
    lfp_index2time(lfp_TrialRec(trial, 1)) - reftime );
reltrialend = min(reltrialend, ...
    lfp_index2time(lfp_TrialRec(trial, 2))  - reftime );
