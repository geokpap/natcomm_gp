function ts = lfp_filtertimestamps(ts, aligns, offsets)
%badtrials = lfp_filtertimestamps(ts)
%   Applies the same logic as the 'windows' option of lfp_findbadtrials to
% identify members of <ts> that do not belong to the analysis window, and
% removes from the returned list.
%INPUTS
% ts: a vector of timestamps in seconds.
% aligns: cell vector of event IDs to use as alignment events.
% offsets: two-column numeric array containing one row for each element of
%   <aligns> where the first column specifies the window start time
%   relative to the alignment event and the second column species the
%   window end time.
%OUTPUTS
% ts: a copy of the <ts> that was submitted, with all timestamps deleted
%   that do not fall within any of the windows specified by <aligns> and
%   <offsets>.
%EXAMPLES
%NOTES
% Does not use the most efficient strategy, but the code is simple.  At the
% cost of much greater complexity, we could sort <ts> and use pointers to
% bracket the range of <ts> that has to be evaluated with regard to each
% trial.

%$Rev: 403 $
%$Date: 2019-08-27 16:52:36 -0400 (Tue, 27 Aug 2019) $
%$Author: dgibson $

global lfp_TrialIndex lfp_Events

isbadTS = true(size(ts));
for trial = lfp_enabledTrials
    for alignidx = 1:numel(aligns)
        trialevts = lfp_Events( ...
            lfp_TrialIndex(trial,1) : lfp_TrialIndex(trial,2), : );
        evtidx = find(ismember(trialevts(:,2), aligns{alignidx}), 1);
        if ~isempty(evtidx)
            evtTS = trialevts(evtidx, 1);
            winTS = evtTS + offsets(alignidx, :);
            isbadTS(ts >= winTS(1) & ts < winTS(2)) = false;
        end
    end
end
ts(isbadTS) = [];

