function [interval, rawtriginfo, badtrigs] = lfp_findCommonTime(trials, ...
    varargin)
%[interval, rawtriginfo] = lfp_findCommonTime(trials)
%   Returns the time interval shared by all trialnums in <trials> expressed
%   in samples relative to lfp_AlignmentRef, plus the <rawtriginfo> table
%   on which <interval> is based.  Note that lfp_findCommonTime differs
%   from  most lfp_lib functions in that <trials> is used literally as
%   supplied, without processing for disabled trials, etc.  If there is no
%   lfp_AlignmentRef in a given trial, an error is raised unless 'norefOK'
%   is specified; in that case, the trial is excluded from the calculation
%   of interval and the corresponding row of rawtriginfo is all zeros.
%   Also, lfp_findCommonTime does not care about lfp_XLimAll, except in the
%   special case of 'multitrig'.
%       <interval> is expressed in units of sample counts, with negative
%   numbers indicating samples before the reference event.  If there are no
%   trials that have lfp_AlignmentRef, <interval> is [NaN NaN]; it is never
%   empty.  It is defined to be a bug in lfp_findCommonTime if the choice
%   of different options results in different quantization-to-samples in
%   <rawtriginfo>.  I.e., lfp_findCommonTime guarantees that <rawtriginfo>
%   will remain the same across all options invocations except for the
%   'continuous' option.
%INPUTS
% trials:  verbatim list of trials.  It is an error for it to be empty.
%   Note that <trials> is not modified in any way by this function, so it
%   should already be filtered according to selection criteria etc, e.g. by
%   calling as lfp_gatherTrialSpikes(lfp_enabledTrials, ...).
%OUTPUTS
% interval: two element vector of times relative to alignment reference,
%   expressed as numbers of samples.
% rawtriginfo: three column array with one row per trigger (alignment
%   reference), containing [startsample endsample refsample] for each.  The
%   triggers are in chronological order, but there may be rows of [0 0 0]
%   intercalated to mark trials that had no reference if 'norefOK'.
% badtrigs: triggers that were eliminated by 'notrunc'.  <badtrigs> is an
%   index into the original list of all triggers; this is not very useful
%   when using 'multitrig', but when using the default single triggering,
%   trials(badtrigs) contains the trial numbers that were skipped. (Note
%   that 'norefOK' does not eliminate any trials; see <rawtriginfo>).
%OPTIONS
%   'continuous' - avoids all attempts at converting times to sample
%       numbers and returns all values as timestamps in seconds.  Cannot be
%       used with 'recseg'. Required when no CSC data are loaded.
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
%       NOTE:  AS OF 5/6/12, IT APPEARS NEITHER THIS WARNING NOR THIS ERROR
%       ARE ACTUALLY IMPLEMENTED:
%       "Issues a warning if there are two successive instances that are
%       within lfp_XLimAll(2) - lfp_XLimAll(1) of each other.  Raises an
%       error if lfp_XLimAll and <window> are both empty.  <rawtriginfo>
%       now contains one row for each trigger rather than one row for each
%       trial."
%   'nobounds' - overrides both 'recseg' and 'continuous' to completely
%       ignore the whole issue of trial or recseg boundaries.
%   'norefOK' -  skips trials that have no alignment event (i.e. treats
%       them as if they were not enabled, but see <rawtriginfo>).
%   'notrunc', xlimspec - simply ignores any triggers whose inclusion would
%       result in truncation of <interval> to something narrower than
%       <xlimspec>; this means that those triggers do not contribute to
%       <rawtriginfo> or <interval>, and are listed in <badtrigs>.  Another
%       way to state this is that all triggers returned have time windows
%       that are entirely contained between the containing trial's start
%       time and end time; if 'recseg' is also used, then each time window
%       is entirely contained within a continuously recorded time segment,
%       but not necessarily within the containing trial.
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

%$Rev: 389 $
%$Date: 2017-06-12 15:07:03 -0400 (Mon, 12 Jun 2017) $
%$Author: dgibson $

rawtriginfo = [];

global lfp_TrialIndex lfp_Events lfp_AlignmentRef lfp_TrialRec ...
    lfp_SamplePeriod lfp_Samples lfp_ActiveFilenums

if isempty(trials)
    error('lfp_findCommonTime:notrials1', ...
        '<trials> is empty.');
end

argnum = 1;
badtrigs = [];
continuousflag = false;
evtbounds = {};
multitrigflag = false;
noboundsflag = false;
norefOKflag = false;
notruncflag = false;
recsegflag = false;
trigfuncflag = false;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'continuous'
            continuousflag = true;
        case 'evtbounds'
            argnum = argnum + 1;
            evtbounds = reshape(varargin{argnum}, 1, []);
            if ~iscell(evtbounds) || ~isequal(size(evtbounds), [1 2])
                error('lfp_findCommonTime:badevtbounds', ...
                    '''evtbounds'' requires a two-element cell value.');
            end
        case 'multitrig'
            multitrigflag = true;
        case 'notrunc'
            notruncflag = true;
            argnum = argnum +1;
            xlimspec = varargin{argnum};
            if length(xlimspec) ~= 2 || ~isnumeric(xlimspec)
                error('lfp_findCommonTime:badxlimspec', ...
                    '<xlimspec> must be a two-element numeric array');
            end
        case 'nobounds'
            noboundsflag = true;
        case 'norefOK'
            norefOKflag = true;
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

if recsegflag && continuousflag
    error('lfp_findCommonTime:conflict', ...
        '''continuous'' and ''recseg'' cannot be used together.');
end

for trialidx = 1:length(trials)
    trial = trials(trialidx);
    startevtidx = lfp_TrialIndex(trial,1);
    endevtidx = lfp_TrialIndex(trial,2);
    if ~isempty(evtbounds)
        startevtidx = find(...
            ismember(...
            lfp_Events(startevtidx:endevtidx, 2), evtbounds{1} )) ...
            + startevtidx - 1; 
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
    % Find <refs> for this trial, the lfp_AlignmentRef absolute
    % timestamp(s). Note that this may be a scalar or a column vector
    % depending on <multitrigflag>, and a floating point number or an
    % integer depending on <continuousflag>.  (Note that <reftime> is
    % always floating-point.)
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
        if norefOKflag
            rawtriginfo = [
                rawtriginfo
                0 0 0 ]; %#ok<AGROW>
        else
            error('lfp_findCommonTime:noref', ...
                'Trial %d has no reference event', trial);
        end
    else
        if ~multitrigflag
            reftime = reftime(1);
        end
        if continuousflag
            refs = reftime;
        else
            refs = lfp_time2index(reftime);
        end
        % <refs> is now a column vector if multitrig, and a scalar if not.
        % It is floating-point if 'continuous', and positive integer(s) if
        % not.
        if noboundsflag
            startendinfo = [1 numel(lfp_Samples{lfp_ActiveFilenums(1)})];
        elseif recsegflag
            startendinfo = [lfp_TrialRec(trial,1) lfp_TrialRec(trial,2)];
        else
            if continuousflag
                startendinfo = reshape( lfp_Events( ...
                    [lfp_TrialIndex(trial,1) lfp_TrialIndex(trial,2)], 1 ), ...
                    1, [] );
            else
                startendinfo = [lfp_TrialIndex(trial,3) lfp_TrialIndex(trial,4)];
            end
        end
        if multitrigflag
            startendinfo = repmat(startendinfo, length(refs), 1);
        end
        % append the current trial's trigger info to <rawtriginfo>:
        rawtriginfo = [
            rawtriginfo
            startendinfo refs ]; %#ok<AGROW>
    end
end

if notruncflag && ~isempty(rawtriginfo)
    % Remove trials from <rawtriginfo> per 'notrunc'
    if continuousflag
        xlimallvals = xlimspec;
    else
        xlimallvals = round(xlimspec/lfp_SamplePeriod);
    end
    istoonarrow = rawtriginfo(:,1) - rawtriginfo(:,3) >  xlimallvals(1) ...
        | rawtriginfo(:,2) - rawtriginfo(:,3) < xlimallvals(2);
    rawtriginfo(istoonarrow, :) = [];
    badtrigs = find(istoonarrow);
end

if isempty(rawtriginfo)
    warning('lfp_findCommonTime:notrials2', 'No trials meet all criteria.');
    interval = [NaN NaN];
else
    triginfo = rawtriginfo(rawtriginfo(:,3) ~= 0, :);
    if isempty(triginfo)
        interval = [NaN NaN];
    else
        pointsbefore = min(triginfo(:,3) - triginfo(:,1));
        pointsafter = min(triginfo(:,2) - triginfo(:,3));
        interval = [ -pointsbefore pointsafter ];
    end
end

