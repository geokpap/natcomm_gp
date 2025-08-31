function [evtmatrix, evtidx] = lfp_getevtmatrix( ...
    triginfo, winpts, evts2avg)
% Collect event data in the style to which lfp_lib wishes it were
% accustomed (see lfp_getSamples).
%INPUTS
% triginfo: as returned by lfp_findCommonTime.
% winpts: A 2-element vector containing the analysis time window in samples
%   relative to the alignment reference sample.  Only events that fall
%   inside <winpts> with respect to each trigger are returned.  If <winpts>
%   is empty, then all events for each triggered trial are returned.
%   Default: empty.
% evts2avg: integer vector of event IDs to gather.  If empty, then all
%   event IDs are gathered.  Default: empty.
%OUTPUTS
% evtmatrix: as returned by lfp_getSamples.  This is a matrix of almost-raw
%   event data in the 2-column format of lfp_Events, except that timestamps
%   for each row of <evtmatrix> are expressed relative to the trigger
%   sample specified in each row of <triginfo(:,3)>.
% evtidx: for each trigger (i.e. row in <triginfo>), evtidx{trigidx}
%   contains the list of indices into lfp_Events that is used to construct
%   <evtmatrix>.  Note that this differs from the format of lfp_getSpikes
%   Rev 325 (the current version as of this writing).
%OPTIONS
% Options are entirely argument-value driven.  If <evts2avg> is not empty,
% then events are gathered as for the 'evtavg'/'evtavg2' options in
% lfp_getSamples.  <winpts> controls which option is executed: if empty,
% then for each trigger, one copy is returned of all the events that belong
% to the same trial as the trigger.  If non-empty, then for each trigger,
% one copy is returned of all the events that fall within <winpts> with
% respect to the trigger.
%NOTES
% For greatest flexibility, we treat triggers as belonging to the last
% trial whose start sample is less than or equal to the trigger sample
% (i.e. we include the following ITI as part of the trial's time).  Any
% triggers that are before the first trial start do not contribute any
% events to <evtmatrix>, <evtidx>.  We rely here on lfp_findCommonTime's
% guarantee that the triggers in <triginfo> are in chronological order.

%$Rev: 394 $
%$Date: 2019-02-14 18:06:15 -0500 (Thu, 14 Feb 2019) $
%$Author: dgibson $

global lfp_Events lfp_TrialIndex

if nargin < 2
    winpts = [];
end
if nargin < 3
    evts2avg = [];
end

evtmatrix = cell(size(triginfo, 1), 1);  % cell array to avoid copying
evtidx = cell(size(triginfo, 1), 1);
reftime = lfp_index2time(triginfo(:,3));
trialnum = 0;
for trigidx = 1:size(triginfo, 1)
    if isempty(winpts)
        % Just take all the events that belong to <trigidx>'s trial.
        % First, skip over all the trials that are too early:
        while trialnum < size(lfp_TrialIndex,1) && ...
                triginfo(trigidx,3) >= lfp_TrialIndex(trialnum+1, 3)
            trialnum = trialnum + 1;
        end
        % If <trialnum> is still 0, that means that either there are no
        % trials (that would be an error), or the current trig sample is
        % less than the first sample in the first trial (and so we ignore
        % the trig).
        if trialnum > 0
            evtidx{trigidx} = (lfp_TrialIndex(trialnum, 1) : ...
                lfp_TrialIndex(trialnum, 2))';
        end
    else
        % potentially slow brute force search for events in <trigidx>'s
        % time window:
        evtidx{trigidx} = reshape( find( lfp_Events(:,1) >= ...
            lfp_index2time(winpts(1) + triginfo(trigidx,3)) ...
            & lfp_Events(:,1) <= ...
            lfp_index2time(winpts(2) + triginfo(trigidx,3)) ), [], 1 );
    end
    if ~isempty(evts2avg)
        evts2del = ~ismember(lfp_Events(evtidx{trigidx}, 2), evts2avg);
        evtidx{trigidx}(evts2del) = [];
    end
    evtidx{trigidx} = reshape(evtidx{trigidx}, [], 1);
    % Save events, with reftime subtracted from timestamps:
    evtmatrix{trigidx} = [
        lfp_Events(evtidx{trigidx}, 1) - reftime(trigidx) ...
        lfp_Events(evtidx{trigidx}, 2)
        ];
end
evtmatrix = cell2mat(evtmatrix); % concatenate all the cells in the array
end

