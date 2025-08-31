function timediff = lfp_evtTimeDiff(trials, evtids)
%timediff = lfp_timediff(trials, evtids)
% For each trial in <trials>, regardless of whether it is enabled or not
% (see lfp_enabledTrials): Using the timestamp of event <evtids(end)> as
% the zero reference, computes the relative time of each other event listed
% in <evtids> and returns the relative times. Each trial in <trials> must
% either have all of the <evtids> or be missing the last one (perhaps along
% with others).  If there is more than one event <evtids(end)>, the first
% one is used.

%$Rev: 301 $
%$Date: 2013-05-03 16:19:03 -0400 (Fri, 03 May 2013) $
%$Author: dgibson $

global lfp_Events lfp_TrialIndex

timediff = NaN(length(trials), length(evtids)-1);
for trialidx = 1:length(trials)
    trial = trials(trialidx);
    trialevents = lfp_Events( ...
        lfp_TrialIndex(trial,1):lfp_TrialIndex(trial,2), : );
    finalevtidx = find(trialevents(:,2) == evtids(end));
    if isempty(finalevtidx)
        continue
    end
    finalevtTS = trialevents(finalevtidx,1);
    for evtnum = 1:length(evtids)-1
        evtidx = find(trialevents(:,2) == evtids(evtnum));
        if isempty(evtidx)
            continue
        end
        timediff(trialidx, evtnum) = trialevents(evtidx(1),1) - ...
            finalevtTS(1);
    end
end
