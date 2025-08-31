function lfp_createEvents(fhandle, eventID, varargin)
%LFP_CREATEEVENTS inserts new events with an arbitrary TTL ID number.
%lfp_createEvents(fhandle, eventID)
%lfp_createEvents(fhandle, eventID, varargin)

%lfp_createEvents(fhandle, eventID)
% <fhandle> is a function handle that points to a function that returns a
% vector of event timestamps in seconds. The timestamps returned are added
% to lfp_Events with <eventID> in the event ID column, in such a way as to
% preserve its timestamp-sorted order.  If multiple events end up having
% the same time stamp, their order is undefined and may change at any time.
% Note that this implies that there is no way of knowing ahead of time
% whether events that have the same timestamp as lfp_NominalTrialStart or
% lfp_NominalTrialEnd will be included in the trial or not (because that
% depends on whether they occur between the trial start and trial end
% events in lfp_Events).
% Expands lfp_EventNames and lfp_EventColors if necessary to accomodate the
% new <eventID>.  It is therefore suggested to avoid using very large
% values for <eventID> (like more than a million). Since Cheetah cannot
% record event IDs greater than 65535, it may be a good idea to use
% <eventID>s greater than 65535; however, then you will be stuck with at
% least 4 MB of ''s in lfp_EventColors and lfp_EventNames. The new entries
% in lfp_EventNames are empty, and the new entries in lfp_EventColors are
% initialized to ''. <eventID> can also be a list of event IDs, in which
% case fhandle must return a cell vector of timestamp vectors, with one
% cell for each event ID.
%lfp_createEvents(fhandle, eventID, varargin)
% Some <fhandles>s may require additional parameters, in which case they
% are put in the place of <varargin>.
%NOTE
% You can insert events at arbitrary literal times by using an anonymous
% function handle, e.g.
%   x=2921.622205; lfp_createEvents(@(a) [a+1 a+2 a+3 a+4], 13, x);
% inserts event 13s at t = 2922.622205, 2923.622205, 2924.622205, and
% 2925.622205.  You can also say
%   lfp_createEvents(@(a) a, sumdumID, []);
% in order to invoke sorting lfp_Events and re-creating lfp_TrialIndex and
% (if necessary) lfp_SelectedTrials without actually adding events,
% where <sumdumID> is just a dummy integer to keep the syntax happy.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

lfp_declareGlobals;

if ~strcmp(class(fhandle), 'function_handle')
    error('lfp_createEvents:badfhandle', ...
        '<fhandle> must be a function handle');
end
if ~strcmp(class(eventID), 'double') ...
        || ~(size(eventID, 1) == 1 || size(eventID, 2) == 1) ...
        || ~(all(fix(eventID) == eventID))
    err(lfp_createEvents:badeventID', ...
        '<eventID> must be an integer vector');
end

disp('Computing event times...');
try
    new_events_array = feval(fhandle, varargin{:});
catch
    rethrow(lasterror);
end

if ~isempty(new_events_array)
    for idx = 1:length(eventID)
        if length(eventID) > 1
            new_events = reshape(new_events_array{idx}, [], 1);
        else
            new_events = reshape(new_events_array, [], 1);
        end
        if (eventID(idx) > length(lfp_EventNames))
            old_EventNames = lfp_EventNames;
            lfp_EventNames = cell(eventID(idx), 1);
            lfp_EventNames(1:length(old_EventNames)) = old_EventNames;
            old_EventColors = lfp_EventColors;
            lfp_EventColors = cell(eventID(idx), 1);
            lfp_EventColors(1:length(old_EventColors)) = old_EventColors;
            lfp_EventColors(max(1,length(old_EventColors)):end) = {''};
            lfp_SelectedEventIDs(end+1:length(lfp_EventColors)) = true;
        end
        % append the new events to lfp_Events:
        lfp_Events = [lfp_Events
            [new_events repmat(eventID(idx), size(new_events))] ];
    end
end

lfp_Events = sortrows(lfp_Events);
disp('Recalculating lfp_TrialIndex...');
lfp_TrialIndex = lfp_createTrialIndex(lfp_Events);
if ~isempty(lfp_SamplePeriod)
    [startSampleIndex, endSampleIndex] = ...
        lfp_addTIcols3n4(lfp_TrialIndex(:,1), lfp_TrialIndex(:,2), lfp_Events);
    lfp_TrialIndex = [lfp_TrialIndex startSampleIndex endSampleIndex];
end
if length(lfp_SelectedTrials) ~= size(lfp_TrialIndex,1)
    % The number of trials has changed, rendering the old
    % lfp_SelectedTrials meaningless and out-of-sync.
    lfp_SelectedTrials = repmat(true, 1, size(lfp_TrialIndex,1));
end

arginstring = '';
if length(varargin) > 0
    arginstring = dg_thing2str(varargin{1});
    for argnum = 2:length(varargin)
        arginstring = [ arginstring ', ' dg_thing2str(varargin{argnum}) ];
    end
end
lfp_log(sprintf('Created new events %s: %s(%s)', mat2str(eventID), ...
    func2str(fhandle), arginstring ));
disp('Done.');
