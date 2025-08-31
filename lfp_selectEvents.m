function result = lfp_selectEvents(trials, varargin)
%result = lfp_selectEvents(trials)
%result = lfp_selectEvents(..., 'evtbounds', {startIDs stopIDs})
%result = lfp_selectEvents(..., 'multitrig')
%result = lfp_selectEvents(..., 'session', sessionname)

% Returns a column vector of the timestamps of events that are selected
% according to the following criteria.  First, the event must be contained
% within one of the specified <trials>, i.e. at or after that trial's
% lfp_NominalTrialStart and at or before the trial's lfp_NominalTrialEnd
% (except when using 'evtboundoffsets'). This implies that a timestamp will
% always be returned if lfp_AlignmentRef is equal to either
% lfp_NominalTrialStart or lfp_NominalTrialEnd. As usual, if any of
% <trials> is not selected in lfp_SelectedTrials, or is listed in
% lfp_BadTrials, that trial is skipped.  If <trials> is [], its value is
% taken to be all trials. <trials> can be specified as trial numbers or as
% Unique Trial IDs.
%
% Second, the event must be an lfp_AlignmentRef event, meaning that
% its event ID must be one of the value listed in lfp_AlignmentRef. 
%
% Third, various options can be invoked (these are copied wherever possible
% from the options to lfp_disp):
%
% OPTIONS
% 'evtboundoffsets', [startoffset stopoffset]
%   Offsets in seconds to add to 'evtbounds' timestamps.  Default = 0.
%   Note that the use of this option can result in timestamps being
%   returned that belong to trials other than the one containing the
%   relevant 'evtbounds' (i.e. if the offsets are large enough to go beyond
%   the trial boundaries).  When that happens, lfp_time2trial is no longer
%   appropriate for determining which trials returned timestamps and which
%   didn't.
% 'evtbounds', {startIDs stopIDs}
%   Sets the range of time to search for alignment events to something
%   narrower than the entire trial, specified in terms of other events.
%   <startIDs> is a list of alternative event IDs to use as the start of the
%   time range, and <endIDs> is a list of alternative event IDs to use as
%   the end of the time range.  The start event is the earliest event in the
%   trial that is a member of <startIDs>; the end event is the earliest
%   event in the trial AFTER the start event that is a member of <endIDs>.
%   It is an error if any trial lacks a start event or end event.
% 'index'
%   Instead of timestamps, returns indices into lfp_Events.
% 'norefOK'
%   Simply skips any trials that do not have a reference event. Default is
%   to raise an error.
% 'multitrig'
%   If there are multiple occurrences of the event in a trial, all
%   occurences are returned.  The default behavior is to return ONLY the
%   FIRST occurrence of lfp_AlignmentRef per trial.
% 'nomultitrig'
%   Raises an error instead of simply using the first instance of the event
%   when multiple instance are found in a trial.
% 'session', sessionname
%   Supplies a default session name to use when specifying Unique Trial IDs
%   instead of internal trial numbers.

%$Rev: 387 $
%$Date: 2016-12-19 16:48:59 -0500 (Mon, 19 Dec 2016) $
%$Author: dgibson $

global lfp_Events
lfp_declareGlobals; % in case mlint missed some

% Process options:
argnum = 1;
evtboundoffsets = [0 0];
evtbounds = {};
indexflag = false;
multitrigflag = false;
nomultitrigflag = false;
norefOKflag = false;
session = '';
while argnum <= length(varargin)
    if ~(ischar(varargin{argnum}) || numel(varargin{argnum}) == 1)
        error('lfp_selectEvents:option', ...
            'Bad option: %s', dg_thing2str(varargin{argnum}));
    end
    switch varargin{argnum}
        case 'evtboundoffsets'
            argnum = argnum + 1;
            evtboundoffsets = varargin{argnum};
        case 'evtbounds'
            argnum = argnum + 1;
            evtbounds = varargin{argnum};
        case 'index'
            indexflag = true;
        case 'multitrig'
            multitrigflag = true;
        case 'nomultitrig'
            nomultitrigflag = true;
        case 'norefOK'
            norefOKflag = true;
        case 'session'
            argnum = argnum + 1;
            session = varargin{argnum};
    end
    argnum = argnum + 1;
end

if nomultitrigflag && multitrigflag
    error('lfp_selectEvents:trig', ...
        '''multitrig'' and ''nomultitrig'' are mutually exclusive; pick one.');
end

% Process <trials> param:
if ischar(trials)
    trials = lfp_parseTrialStr(trials, session);
end
if isempty(trials)
    trials = 1:size(lfp_TrialIndex,1);
end
if any(trials > size(lfp_TrialIndex,1) | trials < 1)
    error('lfp_selectEvents:badtrialnum', ...
        [ 'You have specified trials that do not exist: ' ...
            num2str(trials( ...
            trials > size(lfp_TrialIndex,1) | trials < 1 )) ]);
end
trials = lfp_enabledTrials(trials);

% Select the events and accumulate the timestamps:
result = [];
for trial = trials
    startevtidx = lfp_TrialIndex(trial,1);
    endevtidx = lfp_TrialIndex(trial,2);
    if ~isempty(evtbounds)
        startevtidx = find(...
            ismember(...
            lfp_Events(startevtidx:endevtidx, 2), evtbounds{1} )) ...
            + startevtidx - 1;
        if isempty(startevtidx)
            error('lfp_selectEvents:evtbounds1', ...
                'Trial %d has no ''evtbounds'' start event', ...
                trial );
        else
            startevtidx = startevtidx(1);
        end
        endevtidx = find(...
            ismember(...
            lfp_Events(startevtidx:endevtidx, 2), evtbounds{2} )) ...
            + startevtidx - 1;
        if isempty(endevtidx)
            error('lfp_selectEvents:evtbounds2', ...
                'Trial %d has no ''evtbounds'' end event', ...
                trial );
        else
            endevtidx = endevtidx(1);
        end
        if ~isequal(evtboundoffsets, [0 0])
            startTS = lfp_Events(startevtidx,1) + evtboundoffsets(1);
            endTS = lfp_Events(endevtidx,1) + evtboundoffsets(2);
            if evtboundoffsets(1) < 0
                % add events until we've gone too far
                while startevtidx > 0 && lfp_Events(startevtidx,1) > ...
                        startTS
                    startevtidx = startevtidx - 1;
                end
                % back up one
                startevtidx = startevtidx + 1; 
            end
            if evtboundoffsets(1) > 0
                % remove events until the first event qualifies
                while startevtidx < size(lfp_Events,1) && ...
                        lfp_Events(startevtidx,1) < ...
                        startTS
                    startevtidx = startevtidx + 1;
                end
            end
            if evtboundoffsets(2) < 0
                % remove events until the last event qualifies
                while endevtidx > 1 && lfp_Events(endevtidx,1) > ...
                        endTS
                    endevtidx = endevtidx - 1;
                end
            end
            if evtboundoffsets(2) > 0
                % add events until we've gone too far
                while endevtidx <= size(lfp_Events,1) && ...
                        lfp_Events(endevtidx,1) < ...
                        endTS
                    endevtidx = endevtidx + 1;
                end
                % back up one
                endevtidx = endevtidx - 1;
            end
        end
    end
    trialevents = lfp_Events(startevtidx : endevtidx, :);
    trialresult = [];
    refidx = find(ismember(trialevents(:,2), lfp_AlignmentRef));
    if indexflag
        trialresult(:,1) = refidx + startevtidx - 1;
    else
        trialresult(:,1) = trialevents(refidx, 1);
    end
    if isempty(trialresult)
        if norefOKflag
            continue
        else
            error('lfp_selectEvents:noref2', ...
                'Trial %d has no reference event', ...
                trial );
        end
    elseif nomultitrigflag && length(trialresult) > 1
        error('lfp_selectEvents:multitrig', ...
            'Trial %d contains %d alignment events', ...
            trial, length(trialresult));
    elseif ~multitrigflag
        trialresult = trialresult(1);
    end
    result = [result; trialresult]; %#ok<AGROW>
end
