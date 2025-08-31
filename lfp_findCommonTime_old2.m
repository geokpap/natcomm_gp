function [interval, rawtrialinfo] = lfp_findCommonTime(trials, varargin)
%[interval, rawtrialinfo] = lfp_findCommonTime(trials)
%   Returns the time interval shared by all trialnums in <trials> expressed
%   in samples relative to lfp_AlignmentRef, plus the <rawtrialinfo> table
%   on which <interval> is based, containing [startsample endsample
%   refsample] for each trialnum in <trials>.  If there is no
%   lfp_AlignmentRef, then refpoint for that trial is reported as 0, and
%   the trial is excluded from the calculation of <interval>.  <interval>
%   is expressed in units of sample counts, with negative numbers
%   indicating samples before the reference event.
%OPTIONS
%   'evtbounds', {startIDs endIDs} - Sets the range of time to search for 
%       alignment events to something narrower than the entire trial,
%       specified in terms of other events. <startIDs> is a list of
%       alternative event IDs to use as the start of the time range, and
%       <endIDs> is a list of alternative event IDs to use as the end of
%       the time range.  The start event is the earliest event in the trial
%       that is a member of <startIDs>; the end event is the earliest event
%       in the trial AFTER the start event that is a member of <endIDs>.
%   'multitrig' - Accumulates data over multiple instances of the alignment
%       event per trial.  (The default is to use only the first instance.)
%       Issues a warning if there are two successive instances that are
%       within lfp_XLimAll(2) - lfp_XLimAll(1) of each other.  Raises an
%       error if lfp_XLimAll and <window> are both empty.  <rawtrialinfo>
%       now contains one row for each trigger rather than one row for each
%       trial.
%	'recseg' - The bounds are the enclosing recorded segments rather than 
%       start and end of trial.
%   'trigfunc', funcHandle, args - WARNING: COPIED FROM LFP_SPIKEANALYSIS
%       BUT NOT YET TESTED!!!
%       This option provides a way to customize the selection of events
%       that are used as alignment events. <funcHandle> is a handle to a
%       function that accepts as its first argument a two-element vector of
%       indices into lfp_Events, representing the start and end events
%       between which to find triggers, and returns a column vector of
%       timestamps to use as reference times (i.e. triggers); this vector
%       does not have to be chronologically ordered.  If 'multitrig' is
%       specified, all returned triggers are used; if not, then only the
%       first in the list is used.  Note that if the list of triggers is
%       not chronologically ordered, then the
%       warning('lfp_findCommonTime:fastmultitrig'...) may miss some
%       overlapping time windows.  <args> represents ALL of the rest of the
%       arguments given to lfp_findCommonTime, and consequently this option
%       MUST be specified last in the argument list. <args> are appended to
%       the argument list that is passed to <funcHandle>, making it
%       possible to pass additional data such as eye movement tables to
%       <funcHandle>.

%$Rev: 234 $
%$Date: 2011-06-06 00:18:50 -0400 (Mon, 06 Jun 2011) $
%$Author: dgibson $

lfp_declareGlobals;
argnum = 1;
evtbounds = {};
multitrigflag = false;
recsegflag = false;
trigfuncflag = false;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'evtbounds'
            argnum = argnum + 1;
            evtbounds = varargin{argnum};
        case 'multitrig'
            multitrigflag = true;
        case 'recseg'
            recsegflag = true;
        case 'trigfunc'
            trigfuncflag = true;
            argnum = argnum +1;
            trigfuncH = varargin{argnum};
            trigfuncArgs = varargin(argnum+1:end);
            argnum = length(varargin);
        otherwise
            error('lfp_findCommonTime:badoption', ...
                ['The option "' varargin{argnum} '" is not recognized.'] );
    end
    argnum = argnum + 1;
end

rawtrialinfo = [];

for trialidx = 1:length(trials)
    trial = trials(trialidx);
    % find 'reftime', the lfp_AlignmentRef absolute timestamp (note
    % that this may be a scalar or a column vector):
    startevtidx = lfp_TrialIndex(trial,1);
    endevtidx = lfp_TrialIndex(trial,2);
    if ~isempty(evtbounds)
        startevtidx = find(...
            ismember(...
            lfp_Events(startevtidx:endevtidx, 2), evtbounds{1} )) ...
            + startevtidx - 1; %#ok<NODEF>
        if isempty(startevtidx)
            error('lfp_findCommonTime:evtbounds1', ...
                'Trial %d has no ''evtbounds'' start event', ...
                trial );
        else
            startevtidx = startevtidx(1);
        end
        endevtidx = find( ...
            ismember(lfp_Events(startevtidx:endevtidx, 2), ...
            evtbounds{2} )) ...
            + startevtidx - 1;
        if isempty(endevtidx)
            error('lfp_findCommonTime:evtbounds2', ...
                'Trial %d has no ''evtbounds'' end event', ...
                trial );
        else
            endevtidx = endevtidx(1);
        end
    end
    if trigfuncflag
        reftime = feval(trigfuncH, ...
            [startevtidx endevtidx], trigfuncArgs{:} ); 
    else
        trialevents = lfp_Events(startevtidx : endevtidx, :);
        reftime = trialevents( ...
            ismember(trialevents(:,2), lfp_AlignmentRef), ...
            1 );
    end
    if isempty(reftime)
        refpoint = 0;
    else
        if multitrigflag
            refpoint = lfp_time2index(reftime);
        else
            refpoint = lfp_time2index(reftime(1));
        end
    end
    % <refpoint> is now a column vector if multitrig, and a scalar if not.
    if recsegflag
        startendinfo = [lfp_TrialRec(trial,1) lfp_TrialRec(trial,2)];
    else
        startendinfo = [lfp_TrialIndex(trial,3) lfp_TrialIndex(trial,4)];
    end
    if multitrigflag
        startendinfo = repmat(startendinfo, length(refpoint), 1);
    end
    rawtrialinfo = [
        rawtrialinfo
        startendinfo refpoint ]; %#ok<AGROW>
end

trialinfo = rawtrialinfo(rawtrialinfo(:,3) ~= 0, :);
pointsbefore = min(trialinfo(:,3) - trialinfo(:,1));
pointsafter = min(trialinfo(:,2) - trialinfo(:,3));
interval = [ -pointsbefore pointsafter];

