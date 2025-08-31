function lfp_showEvents(a, b, evtlist, reftime)
%lfp_showEvents(a, b)
% Nicely formats and displays the events between timestamp a and timestamp
% b.  If b is not given or is empty, then the time interval is the entire
% trial whose trial number is round(a), or if a is an array, then every
% trial in round(a).  If evtlist is given and not empty, it is a list of
% event IDs to display; otherwise all events are displayed. If reftime is
% given and not empty, it is subtracted from absolute timestamps.

%$Rev: 326 $
%$Date: 2014-05-29 16:21:14 -0400 (Thu, 29 May 2014) $
%$Author: dgibson $

global lfp_TrialIndex lfp_Events lfp_EventNames

if nargin < 4
    reftime = [];
end
if nargin < 3
    evtlist = [];
end
if nargin < 2 || isempty(b)
    trials = round(a);
    for tridx = 1:numel(trials)
        a(tridx) = lfp_Events(lfp_TrialIndex(trials(tridx),1), 1);
        b(tridx) = lfp_Events(lfp_TrialIndex(trials(tridx),2), 1);
    end
end
if isempty(reftime)
    reftime = 0;
end
for aidx = 1:numel(a)
    if isempty(evtlist)
        eix = find(lfp_Events(:,1) >= a(aidx) & lfp_Events(:,1) <= b(aidx));
    else
        eix = find(lfp_Events(:,1) >= a(aidx) ...
            & lfp_Events(:,1) <= b(aidx) ...
            & ismember(lfp_Events(:,2), evtlist) );
    end
    for k = reshape(eix, 1, [])
        if lfp_Events(k,2) < 1
            eventname = 'bogous';
        else
            eventname = lfp_EventNames{lfp_Events(k,2)};
        end
        disp(sprintf('%.6f %8d %s', lfp_Events(k,1) - reftime, ...
            lfp_Events(k,2), eventname ));
    end
end
    