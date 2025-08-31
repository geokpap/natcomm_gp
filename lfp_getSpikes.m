function [spikes, triginfo, evtidx, evtmatrix, noreftrialidx] = ...
    lfp_getSpikes(trials, clustnums, window, varargin)
%[spikes, triginfo, evtidx, evtmatrix, noreftrialidx] = ...
%    lfp_gatherTrialSpikes(trials, clustnums, win)
% Returns a cell array of data, along with metadata that were used in
% constructing it.  See lfp_getSamples for discrete-time CSC data.
%INPUTS
% clustnum:  cluster to gather from.  It is an error for
%   it to be empty.
% trials:  verbatim list of trials to gather from.  It is an error for it
%   to be empty.  Note that <trials> is not modified in any way by this
%   function, so it should already be filtered according to selection
%   criteria etc, e.g. by calling as
%   lfp_gatherTrialSpikes(lfp_enabledTrials, ...).
% win:  determines <interval> together with lfp_XLimAll.
%OUTPUTS
% spikes: a cell array of column vectors of spike times relative to reftime
%   (which is 0 if ref event does not exist in trial). It is in trials x
%   triggers x clustnums form, where "triggers" refers to multiple
%   occurences of the reference event per trial.  Unused cells are simply
%   left empty, so it may potentially end up containing many empty cells if
%   a small number of trials have many more triggerings than the rest.
%   There is no simple way to determine whether an empty cell represents a
%   trigger that contained no spikes or the absence of a trigger, but
%   size(triginfo,1) is always the actual number of triggers that occurred.
% triginfo: a three-column array with one row for each trigger.  When not
%   using 'multitrig', that means one row for each trial.  The columns
%   contain respectively absolute start time, absolute end time, and
%   absolute trigger time for collection of one trigger's set of spikes.
% evtidx: if 'evts' option is specified, evtidx{trialidx, trigidx} contains
%   the corresponding list of indices into lfp_Events.  Otherwise, <evtidx>
%   is [].
% evtmatrix: contains events in lfp_Events format, but with time relative
%   to the reference event for each trigger.  The number of rows depends on
%   the events that occur in each trial.  It is in overall chronological
%   order within the session, but this means that the relative timestamps
%   in column 1 are NOT in ascending order if there is more than one trial.
% noreftrialidx: a list of unique <trialidx> values for the trials that
%   were skipped by use of 'norefOK'.  Empty if 'norefOK' is not used.  The
%   trialnums of the skipped trials would be trials(noreftrialidx).
%OPTIONS
% 'evtavg', evts2avg - Works like 'evtavg' in lfp_spec.  Collects times of
%   events whose IDs are in <evts2avg> and plots an event marker at the
%   median event time with clickable info.  Due to roundoff errors, results
%   are only good to about 1 part in 1e+4.  Does nothing when <type> is
%   'place' or 'phase'.
% 'evtavg2', evts2avg - Same as 'evtavg', except only events that fall
%   within the trial boundaries AND the display window are used.
% 'evtbounds', evtbounds -  Sets the range of time to search for alignment
%   events to something narrower than the entire trial, specified in terms
%   of other events. <startIDs> is a list of alternative event IDs to use
%   as the start of the time range, and <endIDs> is a list of alternative
%   event IDs to use as the end of the time range.  The start event is the
%   earliest event in the trial that is a member of <startIDs>; the end
%   event is the earliest event in the trial AFTER the start event that is
%   a member of <endIDs>.
% 'evts' - causes collection of data returned as <evtidx>.
% 'meanfr' - for use by 'meanfr' option of lfp_spikeAnalysis.
% 'multitrig' - accumulates data over multiple instances of the alignment
%   event per trial.  (The default is to use only the first instance.)
%   Issues a warning if there are two successive instances that are within
%   lfp_XLimAll(2) - lfp_XLimAll(1) of each other.  Raises an error if
%   lfp_XLimAll and <window> are both empty.
% 'norefOK' - simply skips any trials trials that do not have a reference
%   event.
% 'place' - for use by 'place' option of lfp_spikeAnalysis.
% 'trigfunc', funcHandle, args - this option provides a way to customize
%   the selection of events that are used as alignment events. <funcHandle>
%   is a handle to a function that accepts as its first argument a
%   two-element vector of indices into lfp_Events, representing the start
%   and end events between which to find triggers, and returns a column
%   vector of timestamps to use as reference times (i.e. triggers); this
%   vector does not have to be chronologically ordered.  If 'multitrig' is
%   specified, all returned triggers are used; if not, then only the first
%   in the list is used.  Note that if the list of triggers is not
%   chronologically ordered, then the
%   warning('lfp_spikeAnalysis:fastmultitrig'...) may miss some overlapping
%   time windows.  <args> represents ALL of the rest of the arguments given
%   to lfp_gatherTrialSpikes, and consequently this option MUST be
%   specified last in the argument list.  <args> are appended to the
%   argument list that is passed to <funcHandle>, making it possible to
%   pass additional data such as eye movement tables to <funcHandle>.
%   (Hard to say what happens to 'evtavg'/'evtavg2' if <funcHandle> returns
%   triggers that don't belong to any trial.  -DG 20140521)
%NOTES
% Originally named lfp_gatherTrialSpikes; changed for consistency with
% lfp_getSamples.  There was once also a lfp_getTrialSpikes (now
% lfp_getTrialSpikes.bak) that had parameters preRoll, postRoll which
% should eventually be implemented as options here if they are ever needed.

%$Rev: 389 $
%$Date: 2017-06-12 15:07:03 -0400 (Mon, 12 Jun 2017) $
%$Author: dgibson $

global lfp_Spikes lfp_Events lfp_TrialIndex lfp_AlignmentRef

evtavg2flag = false;
evtbounds = {};
evts2avg = [];
evtsflag = false;
meanfrflag = false;
multitrigflag = false;
norefOKflag = false;
placeflag = false;
trigfuncflag = false;

argnum = 1;
while argnum <= length(varargin)
    if ischar(varargin{argnum})
        switch varargin{argnum}
            case 'evtavg'
                argnum = argnum + 1;
                evts2avg = varargin{argnum};
            case 'evtavg2'
                argnum = argnum + 1;
                evts2avg = varargin{argnum};
                evtavg2flag = true;
            case 'evtbounds'
                argnum = argnum + 1;
                evtbounds = varargin{argnum};
            case 'evts'
                evtsflag = true;
            case 'meanfr'
                meanfrflag = true;
            case 'multitrig'
                multitrigflag = true;
            case 'norefOK'
                norefOKflag = true;
            case 'place'
                placeflag = true;
            case 'trigfunc'
                trigfuncflag = true;
                argnum = argnum +1;
                trigfuncH = varargin{argnum};
                trigfuncArgs = varargin(argnum+1:end);
                argnum = length(varargin);
            otherwise
                error('lfp_gatherTrialSpikes:badoption', ...
                    'The option %s is not recognized.', ...
                    dg_thing2str(varargin{argnum}));
        end
    else
        error( 'lfp_gatherTrialSpikes:badoption2', ...
            'The value %s occurs where an option name was expected', ...
            dg_thing2str(varargin{argnum}) );
    end
    argnum = argnum + 1;
end

if isempty(trials)
    error('lfp_gatherTrialSpikes:notrials', ...
        '<trials> is empty.');
end
if isempty(clustnums)
    error('lfp_gatherTrialSpikes:nofiles', ...
        '<clustnums> is empty.');
end

noreftrialidx = [];
triginfo = [];
spikes = cell(length(trials), 1, length(clustnums));
% <evtidx> is like <spikes>, but indices into lfp_Events, and has no
% clutnum index.
evtidx = cell(length(trials), 1);
evtmatrix = []; % for 'evtavg'

for trialidx = 1:length(trials)
    trial = trials(trialidx);
    startingeventrange = lfp_TrialIndex(trial,1) : ...
        lfp_TrialIndex(trial,2);
    % find 'reftime', the lfp_AlignmentRef absolute timestamp (note
    % that this may be a scalar or a column vector):
    startevtidx = lfp_TrialIndex(trial,1);
    endevtidx = lfp_TrialIndex(trial,2);
    if ~isempty(evtbounds)
        startevtidx = find(...
            ismember(...
            lfp_Events(startevtidx:endevtidx, 2), evtbounds{1} )) ...
            + startevtidx - 1;
        if isempty(startevtidx)
            error('lfp_gatherTrialSpikes:evtbounds1', ...
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
            error('lfp_gatherTrialSpikes:evtbounds2', ...
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
        reftime = 0; %#ok<NASGU>
        if norefOKflag
            noreftrialidx = unique([noreftrialidx trialidx]);
            continue
        else
            error('lfp_gatherTrialSpikes:noref', ...
                'Trial %d has no reference event', ...
                trial );
        end
    elseif ~multitrigflag
        reftime = reftime(1);
    end
    % <reftime> is in seconds, and is now a column vector if
    % multitrig, and a scalar if not.
    
    % collect the spikes for this trial (and events conditionally)
    if isempty(window)
        % timerange is a 2-element row vector (note that <window>
        % being empty implies we are not using multitrig):
        if meanfrflag
            timerange = [ lfp_Events(startevtidx, 1) ...
                lfp_Events(endevtidx, 1) ];
        else
            timerange = [ lfp_Events(lfp_TrialIndex(trial,1), 1) ...
                lfp_Events(lfp_TrialIndex(trial,2), 1) ];
        end
    else
        % timerange has as many rows as reftime, and 2 columns:
        timerange = repmat(reftime, 1, 2) + ...
            repmat(window, size(reftime,1), 1);
        if length(reftime) > 1 && ...
                any(reftime(2:end) - reftime(1:end-1) < ...
                window(2) - window(1) )
            warning('lfp_gatherTrialSpikes:fastmultitrig', ...
                'Trial %d contains multitrig events that will cause spikes to be double-counted', ...
                trial );
        end
    end
    % The whole <spkidx1> and <spkidx2> strategy turned out to be
    % self-defeating when I actually measured the performance, so
    % eliminated it.
    for clustidx = 1:length(clustnums)
        clustnum = clustnums(clustidx);
        for trigidx = 1:size(timerange,1)
            spkindices = find( ...
                lfp_Spikes{clustnum} > timerange(trigidx, 1) ...
                & lfp_Spikes{clustnum} < timerange(trigidx, 2) );
            if ~isempty(spkindices) && (trigidx < size(timerange,1))...
                    && ( lfp_Spikes{clustnum}(spkindices(end)) ...
                    < timerange(trigidx+1,1) )
            end
            % Main spike collection
            if placeflag
                spikes{trialidx, trigidx, clustidx} = ...
                    reshape(lfp_Spikes{clustnum}(spkindices), [], 1);
            else
                spikes{trialidx, trigidx, clustidx} = ...
                    reshape( lfp_Spikes{clustnum}(spkindices) ...
                    - reftime(trigidx), [], 1 );
            end
            if evtsflag
                evtidx{trialidx, trigidx} = find( ...
                    lfp_Events(:,1) > timerange(1) ...
                    & lfp_Events(:,1) < timerange(2) );
                if ~isempty(evts2avg)
                    eventrange = startingeventrange;
                    if evtavg2flag
                        if isempty(window)
                            starttime = lfp_Events( ...
                                lfp_TrialIndex(trial,1), 1 );
                            endtime = lfp_Events( ...
                                lfp_TrialIndex(trial,2), 1 );
                        else
                            starttime = max( ...
                                reftime(trigidx) + window(1), ...
                                lfp_Events( ...
                                lfp_TrialIndex(trial,1), 1 ));
                            endtime = min( ...
                                reftime(trigidx) + window(2), ...
                                lfp_Events( ...
                                lfp_TrialIndex(trial,2), 1 ));
                        end
                        eventrange = eventrange( ...
                            lfp_Events(eventrange,1) >= starttime ...
                            & lfp_Events(eventrange,1) <= endtime );
                    end
                    evtidx{trialidx, trigidx} = eventrange;
                    % Convert <evtidx> into <evtmatrix>
                    evtsinrange = lfp_Events(eventrange,:);
                    evts2include = ismember(evtsinrange(:,2), evts2avg);
                    evtmatrix = [evtmatrix
                        [evtsinrange(evts2include,1)-reftime(trigidx) ...
                        evtsinrange(evts2include,2)]
                        ]; %#ok<AGROW>
                end
            end
        end
    end
    % in 'multitrig' mode, triginfo has a row for each trigger
    % event, so it may have many rows per trial (in which case
    % <reftime> is a column vector).
    triginfo = [ triginfo
        timerange-repmat(reftime,1,2) reftime ]; %#ok<AGROW>
end %for trialidx = 1:length(trials)

% Done gathering spikes over trials


