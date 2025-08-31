function ts = lfp_evtTime(trialnum, evtID)
%ts = lfp_evtTime(trialnum, evtID)
% Returns the timestamp of the first event in trial number <trialnum> that
% is a member of the vector of event IDs <eventID>.

%$Rev: 212 $
%$Date: 2011-02-20 00:26:39 -0500 (Sun, 20 Feb 2011) $
%$Author: dgibson $

global lfp_Events lfp_TrialIndex

trialevents = lfp_Events( ...
    lfp_TrialIndex(trialnum,1):lfp_TrialIndex(trialnum,2), : );
evtidx = find(ismember(trialevents(:,2),  evtID), 1);
if isempty(evtidx)
    warning('lfp_evtTime:nse', ...
        'There is no event %s in trial num %d', ...
        dg_thing2str(evtID), trialnum);
    ts = NaN;
else
    ts = lfp_Events(lfp_TrialIndex(trialnum,1) + evtidx - 1, 1);
end
