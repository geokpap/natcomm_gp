function [sampledata, timepts, evtmatrix, evtidx, badtrials, trials, ...
    filenums] = lfp_getSamples(trials, filenums, window, varargin)
%[sampledata, trials, timepts, evtmatrix, evtidx] = lfp_getSamples(...
%   trials, filenums, window)
% Standardized function for formatting a collection of sample data in the
% style to which lfp_lib wishes it were accustomed.  Time windows are
% specified relative to lfp_AlignmentRef.  The default behavior is to
% return data for as much of the specified window as actually exists in the
% recorded segments (see lfp_RecSegments).  If neither <window> nor
% <lfp_XLimAll> specifies a window, then whatever time period is common to
% all trials is returned. If gaps in recording (see lfp_RecSegments) make
% that impossible, then <window> is truncated as needed to return a
% rectangular array containing the same time window from every trial (but
% see 'notrunc' option).
%INPUTS
% trials:  can be a row vector or a column vector, or it can be empty, in
%   which case all enabled trials are used (equivalent to specifying
%   <trials> as find(lfp_SelectedTrials)). If <trials> is a string, then it
%   is interpreted as containing Unique Trial IDs, i.e. the combination of
%   a session name and a trial number as it was originally numbered in that
%   session; the session name component is optional (see lfp_parseTrialStr
%   and lfp_getTrialNum for further details).  Only trials that are enabled
%   (i.e. lfp_SelectedTrials(trial) is 'true' and <trial> is not a member
%   of lfp_BadTrials) are actually included in the output.
% filenums:  specifies the CSC channels to use.  If <filenums> is empty,
%   then all channels are used.  Only filenums that are enabled (i.e.
%   lfp_SelectedFiles(filenum) is 'true') are actually included in the
%   output.
% window:  The time window of samples to return, relative to each alignment
%   event.  By default, there is one alignment event per trial, and it is
%   the first instance of lfp_AlignmentRef in the trial.  It is an error if
%   there is no lfp_AlignmentRef event in one of the trials.  The options
%   'multitrig', 'notrunc', and 'norefOK' can cause more or fewer alignment
%   events to be used. If <window> is empty, then lfp_XLimAll is used.  If
%   lfp_XLimAll is also empty, the entire trial is used.  If <window> is a
%   1x2 array of double, then the value is used directly.  The time window
%   for a given trial may run past the trial start or end, but it will be
%   truncated if necessary to fit within the trial's recorded segment (see
%   lfp_RecSegments).
%OUTPUT
% sampledata:  sample data in (samples, trigs, filenums) format.  <trigs>
%   refers to the individual instances of alignment events found, which by
%   default is one per trial. In the special case where the 'notrunc'
%   option is specified together with empty values for BOTH <window> AND
%   <lfp_XLimAll>, whole-trial-data mode is invoked, in which <sampledata>
%   is a cell array of trial data in (samples, filenums) form.  The only
%   other options that are valid in whole-trial-data mode are 'logical',
%   'norefOK', 'rmdc', 'rmtrend', and 'verbose'.
% timepts:  row vector containing time of each sample relative to Alignment
%   Reference event, except in whole-trial-data mode, in which case
%   <timepts> is a cell array of the same size as <sampledata>, containing
%   a separate row vector of relative time values for each trial.
% evtmatrix: contains all the events for each enabled trial in lfp_Events
%   format, but with timestamps expressed relative to the reference event
%   for each trigger.  When not using 'multitrig', it is in overall
%   chronological order within the session, but the relative timestamps in
%   column 1 are NOT in ascending order if there is more than one trial.
%   Also, when using 'multitrig', there will be as many copies of each
%   trial's events as there are triggers in that trial, and each copy will
%   have different relative timestamps because they are relative to
%   different triggers. This is primarily intended for use by 'evtavg'
%   display options in downstream code.
% evtidx: for each trigger (i.e. column in <sampledata>), evtidx{trigidx}
%   contains the list of indices into lfp_Events that was used to construct
%   <evtmatrix>.  Note that this differs from the format of lfp_getSpikes
%   Rev 325 (the current version as of this writing).
% badtrials: list of trialnums that were eliminated by lfp_findCommonTime.
% trials: the verbatim list of trials as submitted to lfp_findCommonTime
%   after all processing.
% filenums: the verbatim list of CSC channels that were actually used.
%OPTIONS
% 'evtavg', evts2avg - Collects times relative to each trigger of events
%   that belong to the same trial as the trigger and whose IDs are in
%   the integer vector <evts2avg>.  Returns results in <evtmatrix>,
%   <evtidx>.  Note that single-trial relative event times can only be
%   precise to within half a sample period, but averaging over multiple
%   trials will improve the precision.
% 'evtavg2', evts2avg - Same as 'evtavg', except events that fall
%   within <window> (quantized to samples) relative to each trigger are
%   used instead of events that belong to the same trial as the trigger.
%   Note that this differs from the 'evtavg2' behavior of lfp_getSpikes Rev
%   325 (the current version as of this writing), which excludes any events
%   that occur outside "the trial boundaries" (hard to say what that means
%   given that it *does* have 'trigfunc' option).  Also note that if
%   <window> and <lfp_XLimAll> are both empty, invoking inclusion of "the
%   entire trial", what this means for events is "the entire time period
%   common to all trials, quantized to samples" which may not actually
%   include the start and end events even when there is only a single trial
%   being analyzed.
% 'evtbounds' - see documentation for lfp_findCommonTime
% 'logical' - returns data as 'logical' data type.
% 'margintime', margintime - expands the time window spec by <margintime>
%   on both ends. When <window> and lfp_XLimAll are both empty and we would
%   normally invoke whole-trial-data mode, if <margintime> is non-zero then
%   we add it to both ends of the common time across all trials, determine
%   if our recorded segments can accommodate that window, and then truncate
%   the window if necessary to fit into the recorded segments so as to
%   return a rectangular data array. This is useful for e.g. moving window
%   analyses that reduce a range of time points to a single result point.
%   Hard to say how this "should" interact with 'notrunc', so I'm not
%   saying anything until the issue comes up.
% 'multitrig' - see documentation for lfp_findCommonTime.
% 'norefOK' - skips trials that have no alignment event, i.e. treats them
%   as if they were not enabled, but raises a warning.
% 'notrunc' - similar to 'notrunc' in lfp_findCommonTime, except that no
%   window should be given here; <window> or <lfp_XLimAll> is used
%   automatically.  If <window> and <lfp_XLimAll> are both empty and
%   'notrunc' is also specified, whole-trial-data mode is invoked, which
%   returns the entirety of each trial packaged separately into its own
%   cell, and in this case <sampledata> is a cell array.  This is useful
%   for functions like lfp_CSCraster that want to do something coordinated
%   over a set of trials but without truncating the data.
% 'rmdc' - Removes the DC component from each column of <sampledata>, or
%   from each trial in whole-trial-mode.
% 'rmtrend' - Removes the best fit line from each column of <sampledata>,
%   or from each trial in whole-trial-mode.
% 'rmEP' - Removes the average over trials from each column in
%   <sampledata>; does not work with whole-trial-mode (see
%   "error('lfp_getSamples:wholetrial'" in code). 
% 'session', sessionname - supplies a default session name to use when
%   specifying Unique Trial IDs instead of internal trial numbers.
% 'shuffle' - Gathers data as usual, but the triggerings (i.e. the
%   trials, unless 'multitrig' is
%   used) of the second channel are randomly shuffled s.t. no
%   triggering is in its original position (cf. the 'shuffle' option in
%   lfp_spec)
% 'verbose' - sends clues to stdout for use in case of desperation
%NOTES
%   In the case where 'notrunc' is specified together with empty values for
% BOTH <window> and <lfp_XLimAll>, invoking whole-trial-data mode, the
% following options are hard just to define (let alone implement) and raise
% an error:
%       'rmEP'
%       'shuffle'
% and the following options could be useful, but pose complications and are
% thus not implemented, and raise an error:
%       'evtavg'
%       'evtavg2'
%       'evtbounds'
%       'multitrig'
%       'trigfunc'
% Note that in whole-trial-data mode, <timepts> is also cell-formatted to
% match <sampledata>.

%$Rev: 382 $
%$Date: 2016-08-18 15:52:48 -0400 (Thu, 18 Aug 2016) $
%$Author: dgibson $

global lfp_SelectedTrials lfp_TrialIndex lfp_ActiveFilenums ...
    lfp_SelectedFiles lfp_XLimAll ...
    lfp_AlignmentRef lfp_SamplePeriod

evtmatrix = [];
evtidx = [];
badtrials = [];

argnum = 1;
commontimeopts = {};
evtavg2flag = false;
evts2avg = [];
logicalflag = false;
margintime = 0;
multitrigflag = false;
norefOKflag = false;
notruncflag = false;
rmdcflag = false;
rmtrendflag = false;
rmEPflag = false;
session = '';
shuffleflag = false;
verboseflag = false;
wholetrialmode = false;
while argnum <= length(varargin)
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
            commontimeopts = [commontimeopts ...
                {'evtbounds' varargin{argnum}}]; %#ok<*AGROW>
        case 'logical'
            logicalflag = true;
        case 'margintime'
            argnum = argnum + 1;
            margintime = varargin{argnum};
        case 'multitrig'
            multitrigflag = true;
            commontimeopts = [commontimeopts {'multitrig'}];
        case 'norefOK'
            norefOKflag = true;
            commontimeopts = [commontimeopts {'norefOK'}];
        case 'notrunc'
            notruncflag = true;
            if ~isempty(window)
                commontimeopts = [commontimeopts {'notrunc' window}];
            end
        case 'rmdc'
            rmdcflag = true;
        case 'rmtrend'
            rmtrendflag = true;
        case 'rmEP'
            rmEPflag = true;
        case 'session'
            argnum = argnum + 1;
            session = varargin{argnum};
        case 'shuffle'
            shuffleflag = true;
        case 'verbose'
            verboseflag = true;
        otherwise
            error('lfp_getSamples:badoption', ...
                ['The option "' varargin{argnum} ...
                '" is not recognized.'] );
    end
    argnum = argnum + 1;
end

% Resolve option conflicts
if rmdcflag && rmtrendflag
    warning('lfp_getSamples:optconflict', ...
        '''rmtrend'' also removes dc; ignoring ''rmdc''.');
    rmdcflag = false;
end
if shuffleflag && length(filenums) < 2
    warning('lfp_getSamples:shuffle', ...
        'Shuffling applies only to the second channel, which does not exist.');
    shuffleflag = false;
end

% Standard argument processing
if isempty(trials)
    trials = 1:size(lfp_TrialIndex,1);
elseif ischar(trials)
    trials = lfp_parseTrialStr(trials, session);
end
if any(trials > length(lfp_SelectedTrials))
    warning('lfp_getSamples:trials', ...
        'Ignoring trials %s, which are beyond the last trial.', ...
        dg_canonicalSeries(trials(trials > length(lfp_SelectedTrials))) );
    trials(trials > length(lfp_SelectedTrials)) = [];
end
if ~(isnumeric(trials)) || ~all(fix(trials(:)) == trials(:))
    error('lfp_getSamples:badTrials2', ...
        '<trials> must be an integer vector.');
end
trials = lfp_enabledTrials(trials);
if isempty(trials)
    error('lfp_getSamples:notrials', ...
        'There are no enabled trials; please check lfp_BadTrials and trial selection criteria.');
end
if isempty(filenums)
    filenums = sort(lfp_ActiveFilenums);
elseif any(~ismember(filenums, lfp_ActiveFilenums))
    error('lfp_getSamples:noSuchFilenum', ...
        'You requested non-existent file numbers: %s', ...
        dg_canonicalSeries(filenums( ...
        ~ismember(filenums, lfp_ActiveFilenums) )));
end
filenums = filenums(lfp_SelectedFiles(filenums));
if isempty(window)
    window = lfp_XLimAll;
else
    window = reshape(window, 1, []);
    if ~isequal(size(window), [1 2])
        error('lfp_getSamples:badwindow1', ...
            '<window> must be two-element vector.' );
    end
    commontimeopts = [commontimeopts {'recseg'}];
end

if verboseflag
    fprintf('calling lfp_findCommonTime with commontimeopts=%s\n', ...
        dg_thing2str(commontimeopts));
end

% <rawtriginfo> is in (trigidx, [startsample endsample refsample]) format.
[interval, rawtriginfo, badtrigs] = lfp_findCommonTime( trials, ...
    commontimeopts{:} );

if margintime == 0
    % Not having to cope with 'margintime' makes the logic a bit simpler.
    % All the default processing for <window> is now done, so we can use it
    % to find if we are in <wholetrialmode>:
    expwindow = window;
    wholetrialmode = isempty(window) && notruncflag;
    if any([evts2avg, evtavg2flag, rmEPflag, shuffleflag])
        error('lfp_getSamples:wholetrial', ...
            'The options ''evtavg'', ''evtavg2'', ''rmEP'', and ''shuffle'' cannot be used in whole-trial mode.');
    end
else
    if notruncflag && isempty(window)
        % practically speaking, the behavior is the same as when 'notrunc'
        % has not been specified, except that 
        warning('lfp_getSamples:undef', ...
            'The interaction of ''notrunc'' and ''margintime'' is undefined for empty <window>.');
    end
    % We need to expand <window> by <margintime>, even if <window> is
    % empty (in which case we expand the <interval> returned from
    % lfp_findCommonTime).
    if isempty(window)
        % add margintime to both ends of the actual common trial time:
        expwindow = lfp_SamplePeriod * interval + [-margintime margintime];
        if notruncflag
            warning('lfp_getSamples:notrunc', ...
                'The interaction of ''notrunc'' and ''margintime'' is not defined.');
        end
    else
        expwindow = window + [-margintime margintime];
    end
end

% Make sure that alignment events were found:
if isempty(rawtriginfo)
    msg = sprintf( 'No events were found for lfp_AlignmentRef = %s', ...
        dg_thing2str(lfp_AlignmentRef) );
    if norefOKflag
        warning('lfp_getSamples:notrigs', '%s', msg);
        return
    else
        error('lfp_getSamples:notrigs', '%s', msg);
    end
end
% <badtrials> can take a while to compute, so we only do it if the
% value is used:
if nargout >= 5
    if multitrigflag
        gottrig = false(size(trials));
        for tidx = 1:length(trials)
            gottrig(tidx) = ...
                any( rawtriginfo(:,3) >= lfp_TrialIndex(trials(tidx), 3) ...
                & rawtriginfo(:,3) < lfp_TrialIndex(trials(tidx), 4) );
        end
        badtrials = trials(~gottrig);
    else
        badtrials = trials(badtrigs);
    end
end
% If 'norefOK', delete the placeholders for the triggerless trials from
% <rawtriginfo>.
if any(rawtriginfo(:,3)==0)
    msg = sprintf( 'The reference event is missing from trial(s) %s', ...
        dg_canonicalSeries(trials(rawtriginfo(:,3)==0)) );
    if norefOKflag
        warning('lfp_getSamples:noref', '%s', msg);
        del_trial = rawtriginfo(:,3) == 0;
        trials(del_trial) = [];
        triginfo = rawtriginfo(~del_trial, :);
    else
        error('lfp_getSamples:noref', '%s', msg);
    end
else
    triginfo = rawtriginfo;
end

% <triginfo> is not so raw as <rawtriginfo> now: all triggerless trials are
% gone.  Ditto for <trials>.

if wholetrialmode
    % Note that <trials> and <triginfo> have already been appropriately
    % edited as specified by the options.
    [sampledata, timepts] = getWholeTrialData(trials, filenums, triginfo, ...
        logicalflag, rmdcflag, rmtrendflag, verboseflag);
else
    % NOT wholetrialmode.
    % 
    if isempty(window)
        if notruncflag
        else
            interval = lfp_findCommonTime( trials, commontimeopts{:} );
        end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % dfgjhqoi adlfgjhald habdugladufkhg %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Make sure the data actually exist to satisfy <expwindow>:
    if ~any(cellfun(@(a) isequal(a, 'recseg'), commontimeopts))
        commontimeopts = [commontimeopts {'recseg'}];
    end
    
    if expwindow(1) <= lfp_SamplePeriod * (recseginterval(1) - 1) ...
            || expwindow(2) >= lfp_SamplePeriod * (recseginterval(2) + 1)
        % We have failed to accommodate <expwindow>.
        error('lfp_getSamples:margintime2', ...
            'There is not enough recorded time to accommodate the specified margintime of %d s.', ...
            margintime);
    end
    [sampledata, timepts, evtmatrix, evtidx] = getRectangularArrayData( ...
        filenums, expwindow, ...
        interval, triginfo, evts2avg, evtavg2flag, logicalflag, ...
        rmdcflag, rmtrendflag, rmEPflag, shuffleflag, verboseflag, ...
        varargin{:} );
end
end


function [sampledata, timepts, evtmatrix, evtidx] = getRectangularArrayData( ...
    filenums, expwindow, ...
    interval, triginfo, evts2avg, evtavg2flag, logicalflag, ...
    rmdcflag, rmtrendflag, rmEPflag, shuffleflag, verboseflag, varargin)
%INPUTS
% window: in seconds.  If empty, this specifies a fixed time window whose
%   width is equal to the time period in common across all trials as
%   returned by lfp_findCommonTime.

% Find sample indices

global lfp_Samples lfp_TrialIndex ...
    lfp_SamplePeriod ...
    lfp_Events

sampledata = [];
timepts = []; %#ok<NASGU>
evtmatrix = [];
evtidx = []; %#ok<NASGU>
% We define that if the window is less than half a sample wide, that
% means one point (not zero points):
windurpts = round(diff(expwindow)/lfp_SamplePeriod) + 1;
xlimpoints = round(expwindow/lfp_SamplePeriod);
xlimpointserr = xlimpoints * lfp_SamplePeriod - expwindow;
[v, bestidx] = min(abs(xlimpointserr)); %#ok<ASGLU>
if bestidx==1
    xlimpoints = xlimpoints(1) + [0, windurpts - 1];
else
    xlimpoints = xlimpoints(2) + [1 - windurpts, 0];
end
% Set the interval for analysis to xlimpoints clipped by <interval>:
interval(1) = max(xlimpoints(1), interval(1));
interval(2) = min(xlimpoints(2), interval(2));
timepts = (interval(1) : interval(2)) * lfp_SamplePeriod;
idxrange = interval(1):interval(2);
if isempty(idxrange)
    error('lfp_getSamples:nodata', ...
        ['No samples were selected; note that if lfp_XLimAll is\n' ...
        'empty, <window> is clipped to start and end of trial.'])
end
if length(idxrange) < 4
    warning('lfp_getSamples:fewdata', ...
        ['Fewer than 4 samples were selected; note that if lfp_XLimAll is\n' ...
        'empty, <window> is clipped to start and end of trial.'])
end
indices = (repmat(idxrange, size(triginfo(:,3))) ...
    + repmat(triginfo(:,3), size(idxrange)) )';
if verboseflag
    fprintf('Found indices\n');
end

for fileidx = 1:length(filenums)
    if logicalflag
        sampledata(:,:,fileidx) = logical(reshape( ...
            lfp_Samples{filenums(fileidx)}(indices), ...
            size(indices) ));
    else
        sampledata(:,:,fileidx) = reshape( ...
            lfp_Samples{filenums(fileidx)}(indices), ...
            size(indices) );
    end
end
if rmtrendflag
    for trigidx = 1:size(sampledata,2)
        sampledata(:,trigidx,:) = detrend( ...
            squeeze(sampledata(:,trigidx,:)) );
    end
end
if rmdcflag
    sampledata = sampledata - ...
        repmat(mean(sampledata,1),size(sampledata,1),1);
end
if rmEPflag
    sampledata = sampledata - ...
        repmat(mean(sampledata,2),1,size(sampledata,2));
end
if shuffleflag
    if size(sampledata,2) > 1
        switch size(sampledata,2)
            case 2
                trig_perm = [2 1];
            case 3
                % for a 3-trial shuffle, the bad permutations would be
                % 1 2 3
                % 1 3 2
                % 2 1 3
                % 3 2 1
                % and the only good ones would be
                % 2 3 1
                % 3 1 2
                if rand < .5
                    trig_perm = [2 3 1];
                else
                    trig_perm = [3 1 2];
                end
            otherwise
                % for four or more trials, we loop until happy
                trig_vector = 1:size(sampledata,2);
                trig_perm = randperm(size(sampledata,2));
                while any(trig_perm == trig_vector)
                    % For the first badly permuted trial, find a different permuted
                    % trial to trade with s.t. neither one is in its original place
                    % after the trade.
                    badpermidx = find(trig_perm == trig_vector, 1);
                    if badpermidx < size(sampledata,2) / 2
                        % switch with later perm
                        incr = 1;
                    else
                        % switch with earlier perm
                        incr = -1;
                    end
                    perm2switch = badpermidx + incr;
                    for attempt = 1:3
                        if perm2switch == trig_perm(badpermidx) ...
                                || trig_perm(perm2switch) == badpermidx
                            perm2switch = perm2switch + incr;
                        else
                            break
                        end
                    end
                    if attempt > 2
                        error('lfp_getSamples:shuffle2', ...
                            'This is logically impossible');
                    end
                    permval = trig_perm(badpermidx);
                    trig_perm(badpermidx) = trig_perm(perm2switch);
                    trig_perm(perm2switch) = permval;
                end
        end
        sampledata(:,:,2) = sampledata(:,trig_perm,2);
    else
        error('lfp_getSamples:shuffle', ...
            'You cannot shuffle a single trial!');
    end
end

% Having thoroughly resolved everything, we may also collect events.
% Someday there will be 'trigfunc' option that will break the guarantee
% that all trigger events belong to trial time (as opposed to ITI), so we
% treat triggers as belonging to the last trial whose start sample is less
% than or equal to the trigger sample (i.e. we include the following ITI as
% part of the trial's time).  Any triggers that are before the first trial
% start do not contribute any events to <evtmatrix>, <evtidx> when using
% 'evtavg'.  We rely here on lfp_findCommonTime's guarantee that the
% triggers in triginfo are in chronological order.
evtidx = cell(size(triginfo, 1), 1);  % cell array to avoid copying
reftime = lfp_index2time(triginfo(:,3));
trialnum = 0;
for trigidx = 1:size(triginfo, 1)
    if evtavg2flag
        % potentially slow brute force search for events in <trigidx>'s
        % time window:
        evtidx{trigidx} = find( lfp_Events(:,1) >= ...
            lfp_index2time(interval(1) + triginfo(trigidx,3)) ...
            & lfp_Events(:,1) <= ...
            lfp_index2time(interval(2) + triginfo(trigidx,3)) );
    else
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
    end
    if ~isempty(evts2avg)
        % foobar: not attempting this just yet, probably needs revisions.
        warning('lfp_getSamples:evts2avg', ...
            'evtavg and evtavg2 not yet implemented' );
        evtidx{trigidx} = evtidx{trigidx}( ...
            ismember(lfp_Events(evtidx{trigidx}, 2), evts2avg) );
    end
    % Append events specified by <evtidx> onto <evtmatrix>, with reftime
    % subtrected from timestamps:
    evtmatrix = [evtmatrix
        [ lfp_Events(evtidx{trigidx}, 1) - reftime(trigidx) ...
        lfp_Events(evtidx{trigidx}, 2) ]
        ];
end
evtidx = cell2mat(evtidx);
end


function [sampledata, timepts] = getWholeTrialData( ...
    trials, filenums, triginfo, ...
    logicalflag, rmdcflag, rmtrendflag, verboseflag)
% <trials> and <triginfo> are assumed to be in final form here, i.e. any
% trials that have been excluded should be gone from <trials> and
% <triginfo>. 
global lfp_SamplePeriod lfp_TrialIndex lfp_Samples
sampledata = cell(size(trials));
timepts = cell(size(trials));
if verboseflag
    fprintf('lfp_getSamples is running in whole-trial-data mode.\n');
end
for trialidx = 1:length(trials)
    trial = trials(trialidx);
    timepts{trialidx} = (( lfp_TrialIndex(trial, 3) : ...
        lfp_TrialIndex(trial, 4) ) - triginfo(trialidx, 3)) * lfp_SamplePeriod;
    for fileidx = 1:length(filenums)
        trialsamples = reshape(lfp_Samples{filenums(fileidx)}( ...
            lfp_TrialIndex(trial, 3) : lfp_TrialIndex(trial, 4) ), [], 1);
        if logicalflag
            sampledata{trialidx}(:,fileidx) = logical(trialsamples);
        else
            sampledata{trialidx}(:,fileidx) = trialsamples;
        end
    end
    if rmtrendflag
        sampledata{trialidx} = detrend(sampledata{trialidx});
    end
    if rmdcflag
        sampledata{trialidx} = sampledata{trialidx} - ...
            repmat( mean(sampledata{trialidx},1), ...
            size(sampledata{trialidx},1), 1 );
    end
end
end
