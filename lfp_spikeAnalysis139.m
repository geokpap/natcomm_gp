function [result, spikes, ntrigs, hF] = lfp_spikeAnalysis( ...
    type, trials, clustnums, window, varargin)
%LFP_SPIKEANALYSIS provides histograms (PETH and phase) and rasters.
%lfp_spikeAnalysis(type, trials, clustnums, window)
%lfp_spikeAnalysis(type, trials, clustnums)
%lfp_spikeAnalysis(type, trials)
%lfp_spikeAnalysis(type)
%lfp_spikeAnalysis(..., 'bigevts', bigevtcodes)
%lfp_spikeAnalysis(..., 'evtavg', evts2avg)
%lfp_spikeAnalysis(..., 'evtbounds', {startIDs stopIDs})
%lfp_spikeAnalysis(..., figbase)
%lfp_spikeAnalysis(..., 'maxskew', percentage, trialfrac)
%lfp_spikeAnalysis(..., 'multichannel')
%lfp_spikeAnalysis(..., 'session', sessionname)
%lfp_spikeAnalysis(..., 'showtrialnums')
%lfp_spikeAnalysis(..., 'trigfunc')
%lfp_spikeAnalysis(..., 'returnfig')
%lfp_spikeAnalysis(..., 'SeqTypeTitle')
%lfp_spikeAnalysis(..., 'subplots', numrowscols, plotnum)
%lfp_spikeAnalysis(..., 'title', titlestr)
%lfp_ras(...)
%lfp_his(...)
%lfp_his(..., 'aligns', n)
%lfp_his(..., 'binwidth', n)
%lfp_his(..., 'hz')
%lfp_his(..., 'multitrig')
%lfp_his(..., 'norefOK')
%lfp_his(..., 'nbins', n)
%lfp_his(..., 'newbins', n)
%lfp_his(..., 'newbinwidth', n)
%lfp_his(..., 'offset')
%lfp_phase(..., wavefilenum)

%[result, spikes, ntrigs] = lfp_spikeAnalysis( ...
%   type, trials, clustnums, window)
%  For each cluster in <clustnums>, displays all selected trials in
%  <trials>, limited to the time interval <window> around the reference
%  event (lfp_AlignmentRef). If <result> is used, then no figure is
%  displayed and the figure data are returned.
%  <spikes> is a cell array of vectors of spike times relative to
%  reftime (which is 0 if ref event does not exist in trial), trimmed to
%  fit <window> and in case of <type> = 'his', further trimmed to fit
%  actual outer bin boundaries.  <ntrigs> is the number of trigger events
%  (aka "alignment events") that went into the analysis; when 'multitrig'
%  is not used, this is just the same as the number trials.
%
%  <type> must be one of:
%   'cdf' - added 2005 by Joey Feingold; same as 'his', except plots result
%       as a cdf instead of a histogram. Must test to make sure subplots
%       are handled correctly just after around line 310.  <result> is
%       whatever is returned by cdfplot.
%   'his' - invoked by lfp_his(...); <result> is a column vector.  (Note
%       that the actual time interval covered by the histogram turns out to
%       be somewhat less than <window> due to the vagaries of the bin
%       boundary calculations that tend to leave half a binwidth of time
%       uncounted at the beginning and end of <window>.)
%   'meanfr' - must be called with 'evtbounds' option as it calculates
%       the mean firing rate, in hz, between those two events. result{1} is
%       the total spikes/total time and result{2} is a vector containing 
%       the mean firing rate per trial.
%   'phase' - invoked by lfp_phase(...); <result> is a column vector.
%   'place' - same as 'ras', except it searches lfp_FileNames for 'trk X'
%       and 'trk Y' and plots the spike data in video tracker X and Y
%       coordinates instead of in trial-by-trial rasters.
%   'ras' - invoked by lfp_ras(...); result is [], so calling result =
%       lfp_spikeAnalysis('ras', ...) is a CPU-intensive null operation.
%       Trials are plotted in the order submitted, so <trials> = 40:-1:1
%       produces a mirror-symmetric plot to <trials> = 1:40.
%  Any other argument that is not supplied or is specified as [] takes a
%  default value:
%    <window> defaults to lfp_XLimAll; if this value is also [], then
%    each whole trial is used
%    <clustnums> defaults to all clusters
%    <trials> defaults to all trials (except those in lfp_BadTrials)
%  WARNING: the y-axis labels are fragile when <type> = 'ras', so it is
%  recommended that you avoid using the magnifying glass or the ylim
%  function on raster displays.  xlim(...) can be used with abandon.
%  If <type> is 'his', then the following additional options are available:
%    'nbins': must be followed by a number specifying the number of bins
%       in the histogram (default = 100); binwidth is calculated as the
%       maximum integral number of ms such that the time window will
%       accomodate nbins.  If the time window is less than nbins ms, an
%       error is raised. If binwidth is specified, together with nbins, it
%       overrides nbins (the time window does NOT change).
%    'binwidth': must be followed by a number specifying the width in ms of
%       the bins in the histogram (default determined by nbins); does
%       something undesirable with phase histograms, so use 'nbins'
%       instead (the undesirable thing it does is to calculate the number
%       of bins exactly the same way it would for a time histogram)
%    'newbins': same as 'nbins', except that all bin boundary calculations
%       are done in integers until the last step before calling dg_hist at
%       which point they are divided by 1000 (1/1000 is a horribly long
%       infinitely repeating binary number).  Also, all relative spike
%       times are rounded to the "nearest microsecond" before calling
%       dg_hist (the scare quotes are to remind you that 1e-6 is also a
%       horrible infinite repeat binary number, so attempting to round this
%       way inevitably results in truncation error, but at least it should
%       wipe out any accumulated rounding errors in the spike times).  Last
%       of all, 'newbins' is not required to be greater than 3, but only
%       greater than 0.  'offset' has no effect when using 'newbins' or
%       'newbinwidth'.
%           The purpose of all this is to make the binning parameters 
%       easily controllable, whereas roundoff errors in the old binning can
%       sometimes cause bins to be dropped or have a width that's off by 1
%       ms from the intended one.
%	 'newbinwidth': same as 'binwidth', except binning is done as for
%       'newbins'.
%    'hz': counts are divided by binwidth in seconds to yield average
%       firing rate
%    'offset': shifts the bins half a binwidth to the left so that 0 is at
%       a junction between bins rather than in the center of a bin; note
%       that this may cause the first bin to include time where there are
%       no data, and therefore to contain an artifactually low spike count.
%   If <type> is 'phase', then the 'nbins' option works as for 'his'.
%lfp_spikeAnalysis(..., 'aligns', alignmentlist)
%  Creates up to six subplots for up to six different alignment events in
%  one figure.  Each of the event IDs in the cell array <alignmentlist> is
%  used in turn as the alignment event; this is done by actually assigning
%  the value to lfp_AlignmentRef, so any previous value of lfp_AlignmentRef
%  is lost.
%lfp_spikeAnalysis(..., 'bigevts', bigevtcodes)
%  <bigevtcodes> is a cell array whose first column contains column vectors
%  of event IDs whose values are higher than the length of lfp_EventColors,
%  and the second column contains color specs with which to display the
%  events IDs listed in the first column.  The color spec can be a
%  character or a numeric triple.  If an event ID occurs in more than one
%  row of <bigevtcodes>, the first row that matches takes precedence.
%lfp_spikeAnalysis(..., 'evtavg', evts2avg)
%  Works like 'evtavg' in lfp_spec.  Collects times of events whose
%  IDs are in <evts2avg> and plots an event marker at the median event time
%  with clickable info.  Due to roundoff errors, results are only good to
%  about 1 part in 1e+4.  Does nothing when <type> is 'place' or 'phase'.
%lfp_spikeAnalysis(..., 'evtavg2', evts2avg)
%  Same as 'evtavg', except only events that fall within the trial
%  boundaries AND the display window are used.
%lfp_spikeAnalysis(..., 'evtbounds', {startIDs endIDs})
%  Sets the range of time to search for alignment events to something
%  narrower than the entire trial, specified in terms of other events.
%  <startIDs> is a list of alternative event IDs to use as the start of the
%  time range, and <endIDs> is a list of alternative event IDs to use as
%  the end of the time range.  The start event is the earliest event in the
%  trial that is a member of <startIDs>; the end event is the earliest
%  event in the trial AFTER the start event that is a member of <endIDs>.
%lfp_spikeAnalysis(..., figbase)
%  The figure window number is (figbase + clustnum).
%lfp_spikeAnalysis(..., 'maxskew', percentage, trialfrac)
%  Returns the scalar NaN if the numbers of spikes accumulated from each
%  trial are grossly skewed, as defined by <percentage> and <trialfrac>: if
%  the most spike-heavy <trialfrac> fraction of the trials contain more
%  than <percentage> percent of the spikes, then NaN is returned.
%  <percentage> is a number in the interval [0, 100], whereas <trialfrac>
%  is in [0, 1].  Has no effect on plotted figures. Note that NaN will
%  ALWAYS be returned if there are too few spikes (see code).
%  (WARNING: this looks as if it will unceremoniously crash - or worse yet,
%  silently malfunction - if used together with multitrig.  -DG 8/8/07)
%lfp_spikeAnalysis(..., 'multichannel')
%  Plots all units in a single figure window.
%lfp_spikeAnalysis(..., 'norefOK')
%  Simply skips any trials trials that do not have a reference event.
%lfp_spikeAnalysis(..., 'session', sessionname)
%  Supplies a default session name to use when specifying Unique Trial IDs
%  instead of internal trial numbers.
%lfp_spikeAnalysis(..., 'showtrialnums')
%  Labels figure with trial numbers instead of selection rule
%lfp_his(..., 'multitrig')
%  Accumulates data over multiple instances of the alignment event per
%  trial.  (The default is to use only the first instance.)  Issues a
%  warning if there are two successive instances that are within
%  lfp_XLimAll(2) - lfp_XLimAll(1) of each other.  Raises an error if
%  lfp_XLimAll and <window> are both empty.
%lfp_spikeAnalysis(..., 'trigfunc', funcHandle, args)
%  This option provides a way to customize the selection of events that are
%  used as alignment events. <funcHandle> is a handle to a function that
%  accepts as its first argument a two-element vector of indices into
%  lfp_Events, representing the start and end events between which to find
%  triggers, and returns a column vector of timestamps to use as reference
%  times (i.e. triggers); this vector does not have to be chronologically
%  ordered.  If 'multitrig' is specified, all returned triggers are used;
%  if not, then only the first in the list is used.  Note that if the list
%  of triggers is not chronologically ordered, then the
%  warning('lfp_spikeAnalysis:fastmultitrig'...) may miss some overlapping
%  time windows.  <args> represents ALL of the rest of the arguments given
%  to lfp_spikeAnalysis, and consequently this option MUST be specified
%  last in the argument list. <args> are appended to the argument list that
%  is passed to <funcHandle>, making it possible to pass additional data
%  such as eye movement tables to <funcHandle>.
%lfp_spikeAnalysis(..., 'returnfig')
%  Forces creation of a figure when return values are used; if not
%  specified, then hF = [].
%lfp_spikeAnalysis(..., 'SeqTypeTitle')
%  Adds spatial and temporal sequence parameters to the fig title for his
%  and ras plots.
%lfp_spikeAnalysis(..., 'subplots', numrowscols, plotnum)
%  General-purpose subplot feature. The current plot is drawn in the
%  subplot numbered <plotnum> of <numrowscols> suplots in the current
%  figure.
%lfp_spikeAnalysis(..., 'title', titlestr)
%  <titlestr> is a string which is used to title the plot or subplot in
%  place of the normally generated title.

%$Rev: 139 $
%$Date: 2010-07-13 20:32:20 -0400 (Tue, 13 Jul 2010) $
%$Author: dgibson $

lfp_declareGlobals;

result = []; %#ok<NASGU>
spikes = [];
ntrigs = [];
hF = [];

if nargin <4 || isempty(window)
    window = lfp_XLimAll;
elseif ~(strcmp(class(window), 'double') ...
        && isequal(size(window), [1 2]) )
    error('lfp_spikeAnalysis:badwindow', ...
        '<window> must be 1x2 number array.' );
end
if nargin < 3 || isempty(clustnums)
    clustnums = 1:length(lfp_Spikes); %#ok<USENS>
end
if isempty(clustnums)
    error('lfp_spikeAnalysis:noclust', ...
        'No clusters were selected for analysis' );
end

if nargin < 2 || isempty(trials)
    trials = 1:size(lfp_TrialIndex,1);
end

result = [];
rasflag = false;
hisflag = false;
cdfflag = false;
phaseflag = false;
placeflag = false;
subplotsflag = false;
meanfrflag = false;
switch type
    case 'cdf'
        cdfflag = true;
    case 'his'
        hisflag = true;
    case 'phase'
        phaseflag = true;
    case 'place'
        placeflag = true;
    case 'ras'
        rasflag = true;
    case 'meanfr'
        meanfrflag = true;
    otherwise
        error('lfp_spikeAnalysis:badMode', ...
            '<type> must be one of: ras, his.' );
end

if phaseflag
    if ~( (~isempty(varargin)) ...
            && strcmp(class(varargin{1}), 'double') ...
            && (fix(varargin{1}) == varargin{1}) )
        error('lfp_spikeAnalysis:needfilenum', ...
            'You must specify a wave filenum after <window> to use ''phase''' );
    end
    wavefilenum = varargin{1};
    varargin(1) = [];
end

if placeflag
    xchan = find(ismember(lfp_FileNames, 'trk X')); %#ok<USENS>
    ychan = find(ismember(lfp_FileNames, 'trk Y'));
    if isempty(xchan) || isempty(ychan)
        error('lfp_spikeAnalysis:needvideo', ...
            'There must be a ''trk X'' channel and a ''trk Y'' channel in lfp_FileNames');
    end
end

session = '';
trialstyle = lfp_TrialStyle;
percentage = 0; % If 'maxskew', then this will be nonzero
argnum = 1;
bigevts = {};
evtavg2flag = false;
evtbounds = {};
evts2avg = [];
figflag = false;
norefOKflag = false;
offsetflag = false;
alignflag = false;
multichannelflag = false;
multitrigflag = false;
newbinsflag = false;
newbinwidthflag = false;
returnfigflag = false;
SeqTypeTitleflag = false;
nbins = 100;
rawbinwidth = 0;
alignevents = {}; %TF, DG
hzflag = false; %TF 2/21/05
trigfuncflag = false;
titleflag = false;
while argnum <= length(varargin)
    if strcmp(class(varargin{argnum}), 'double')
        if figflag
            error('lfp_spikeAnalysis:multifigbase', ...
                'You have specified figbase more than once.');
        else
            figflag = true;
            figbase = varargin{argnum};
        end
    else
        switch varargin{argnum}
            case 'bigevts'
                argnum = argnum + 1;
                bigevts = varargin{argnum};
                try
                    bigevtIDs = cell2mat(bigevts(:,1));
                catch %#ok<CTCH>
                    error('lfp_spikeAnalysis:badbigevts', ...
                        'Value for <bigevts> is badly formatted.' );
                end
            case 'binwidth'
                argnum = argnum + 1;
                if (argnum > length(varargin)) ...
                        || (~strcmp(class(varargin{argnum}), 'double'))
                    error('lfp_spikeAnalysis:badbinwidth', ...
                        'You must specify a numeric value for binwidth' );
                end
                rawbinwidth = 1.0e-3 * varargin{argnum}; % convert ms to s
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
            case 'hz'
                hzflag = true; %TF 2/21/05
            case 'maxskew'
                if ~isequal(class(varargin{argnum+1}), 'double') ...
                        || ~isequal(class(varargin{argnum+2}), 'double')
                    error('lfp_spikeAnalysis:badmaxskew', ...
                        '''maxskew'' option requires two numeric arguments' );
                end
                argnum = argnum +1;
                percentage = varargin{argnum};
                argnum = argnum +1;
                trialfrac = varargin{argnum};
            case 'multichannel'
                multichannelflag = true;
            case 'multitrig'
                multitrigflag = true;
            case 'nbins'
                argnum = argnum + 1;
                if (argnum > length(varargin)) ...
                        || (~strcmp(class(varargin{argnum}), 'double'))
                    error('lfp_spikeAnalysis:badnbins', ...
                        'You must specify a numeric value for nbins' );
                end
                nbins = varargin{argnum};
                if nbins < 3
                    error('''nbins'' must be at least 3');
                end
            case 'newbins'
                newbinsflag = true;
                argnum = argnum + 1;
                nbins = round(varargin{argnum});
                if nbins < 1
                    error('''newbins'' must be at least 1');
                end
            case 'newbinwidth'
                newbinwidthflag = true;
                argnum = argnum + 1;
                rawbinwidth = round(varargin{argnum});
            case 'norefOK'
                norefOKflag = true;
            case 'offset'
                offsetflag = true;
            case 'aligns' %TF added, DG modified
                argnum = argnum +1;
                alignflag = true;
                alignevents = varargin{argnum};
                oldAlignmentRef = lfp_AlignmentRef; %#ok<NODEF>
            case 'session'
                argnum = argnum + 1;
                session = varargin{argnum};
            case 'showtrialnums'
                trialstyle = 'trialnums';
            case 'trigfunc'
                trigfuncflag = true;
                argnum = argnum +1;
                trigfuncH = varargin{argnum};
                trigfuncArgs = varargin(argnum+1:end);
                argnum = length(varargin);
            case 'returnfig'
                returnfigflag = true;                
            case 'SeqTypeTitle'
                SeqTypeTitleflag = true;                
            case 'subplots'
                subplotsflag = true;
                argnum = argnum + 1;
                numrowscols = varargin{argnum};
                argnum = argnum + 1;
                plotnum = varargin{argnum};                
            case 'title'
                titleflag = true;
                argnum = argnum + 1;
                titlestr = varargin{argnum};
            otherwise
                error('lfp_spikeAnalysis:badoption', ...
                    ['The option "' varargin{argnum} '" is not recognized.'] );
        end
    end
    argnum = argnum + 1;
end

if meanfrflag
    if isempty(evtbounds)
        error('lfp_spikeAnalysis:betweenevtbounds', ...
            '''meanfr'' type must be called with ''evtbounds''');
    end
    if  ~isempty(window)
        error('lfp_spikeAnalysis:betweenevtwindow', ...
            '''meanfr'' type must be called with empty ''window''');
    end
    if nargout == 0
        error('''meanfr'' only returns result and not a figure, you must specify output arg')
    end
    if multitrigflag
        error('''meanfr'' does not work with ''multitrig'' option')
    end
    if length(clustnums) > 1
        warning('lfp_spikeAnalysis:mfrclust', ...
            'result of ''meanfr'' will only reflect the *last* cluster')
    end
end

if nargout > 0 && ~returnfigflag && alignflag
    error('''aligns'' option can only be used on figure displays');
end

if ~strcmp(class(nbins), 'double')
    error('lfp_spikeAnalysis:badnbins', ...
        'Value of nbins must be numeric.');
end

if ischar(trials)
    trials = lfp_parseTrialStr(trials, session);
end
% Apply mask to trials:
trials = lfp_enabledTrials(trials);

if size(trials,1) > 1
    if size(trials,2) > 1
        error('lfp_spikeAnalysis:badTrials1', '<trials> must be an integer vector.');
    end
    trials = trials';
end
% <trials> is now a row vector.
if any(trials > size(lfp_TrialIndex,1) | trials < 1)
    error('lfp_spikeAnalysis:badtrialnum', ...
        [ 'You have specified trials that do not exist: ' ...
            num2str(trials( ...
            trials > size(lfp_TrialIndex,1) | trials < 1 )) ]);
end
if ~(strcmp(class(trials), 'double')) ...
        || ~all(fix(trials) == trials)
    error('lfp_spikeAnalysis:badTrials2', '<trials> must be an integer vector.');
end

if isempty(trials)
    error('lfp_spikeAnalysis:nodata', ...
        'There are no trials selected for processing (check lfp_BadTrials, lfp_SelectedTrials)');
end

if newbinsflag && newbinwidthflag
    error('lfp_spikeAnalysis:badopt4', ...
        'Specify only one of ''newbins'' or ''newbinwidth''.');
end

if multitrigflag
    if ~hisflag
        error('lfp_spikeAnalysis:badopt1', ...
            '''multitrig'' only works on histograms');
    end
    if isempty(window)
        error('lfp_spikeAnalysis:badopt2', ...
            '''multitrig'' requires a pre-set time window');
    end
end

% if don't input any events then set it to the global variable
if ~alignflag
    alignevents{1} = lfp_AlignmentRef;  % DG modified
end

% Array used to return to partly-finished multiplots when using 'aligns':
fighandles(clustnums) = 0;

% run the whole procedure for each alignment ref
for i = 1:length(alignevents)
    if length(alignevents) > 1 || multichannelflag
        multiplotflag = true;
    else
        multiplotflag = false;
    end
    lfp_AlignmentRef = alignevents{i};% DG modified
    
    for clustidx = 1:length(clustnums)
        trialinfo = [];
        norefidx = [];
        cluster = clustnums(clustidx);
        % <spikes> is a cell array of column vectors of spike times
        % relative to reftime (which is 0 if ref event does not exist in
        % trial). It is in trials x triggers form, where "triggers" refers
        % to multiple occurences of the reference event per trial.  Unused
        % cells are simply left empty, so it may potentially end up
        % containing many empty cells if a small number of trials have many
        % more triggerings than the rest.
        spikes = cell(length(trials), 1);
        % <evtidx> is like <spikes>, but indices into lfp_Events.
        evtidx = cell(length(trials), 1);
        spkidx1 = 1;
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
                    + startevtidx - 1; %#ok<NODEF>
                if isempty(startevtidx)
                    error('lfp_spikeAnalysis:evtbounds1', ...
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
                    error('lfp_spikeAnalysis:evtbounds2', ...
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
                reftime = 0;
                if norefOKflag
                    norefidx = unique([norefidx trialidx]);
                    continue
                else
                    if phaseflag || hisflag
                        error('lfp_spikeAnalysis:noref', ...
                            'Trial %d has no reference event', ...
                            trial );
                    else
                        warning('lfp_spikeAnalysis:noref2', ...
                            'Trial %d has no reference event', ...
                            trial );
                    end
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
                    warning('lfp_spikeAnalysis:fastmultitrig', ...
                        'Trial %d contains multitrig events that will cause spikes to be double-counted', ...
                        trial );
                end
            end
            % Find last spike in trial to shorten search time
            if isempty(window)
                if meanfrflag
                    endOfSearchInterval = lfp_Events(endevtidx,1);
                else
                    endOfSearchInterval= lfp_Events(lfp_TrialIndex(trial,2), 1);
                end
            else
                endOfSearchInterval = lfp_Events(endevtidx,1) + window(2);
            end
            spkidx2 = find( lfp_Spikes{cluster}(spkidx1:end) ...
                < endOfSearchInterval ) + spkidx1 - 1;
            if ~isempty(spkidx2)
                spkidx2 = spkidx2(end);
            end
            for trigidx = 1:size(timerange,1)
                spkindices = find( ...
                    lfp_Spikes{cluster}(spkidx1:spkidx2) ...
                    > timerange(trigidx, 1) ...
                    & lfp_Spikes{cluster}(spkidx1:spkidx2) ...
                    < timerange(trigidx, 2) ) ...
                    + spkidx1 - 1;
                % update spkidx1 only if none of current window's spikes
                % will be in next window:
                if ~isempty(spkindices) && (trigidx < size(timerange,1))...
                        && ( lfp_Spikes{cluster}(spkindices(end)) ...
                        < timerange(trigidx+1,1) )
                    spkidx1 = spkindices(end) + 1;
                end
                % Main spike collection
                if placeflag
                    spikes{trialidx, trigidx} = ...
                        reshape(lfp_Spikes{cluster}(spkindices), [], 1);
                else
                    spikes{trialidx, trigidx} = ...
                        reshape( lfp_Spikes{cluster}(spkindices) ...
                        - reftime(trigidx), [], 1 );
                end
                if rasflag || placeflag || ~isempty(evts2avg)
                    evtidx{trialidx, trigidx} = find( ...
                        lfp_Events(:,1) > timerange(1) ...
                        & lfp_Events(:,1) < timerange(2) );
                    if evts2avg
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
            % in 'multitrig' mode, trialinfo has a row for each trigger
            % event, so it may have many rows per trial (in which case
            % <reftime> is a column vector).
            trialinfo = [ trialinfo
                timerange-repmat(reftime,1,2) reftime ]; %#ok<AGROW>
        end %for trialidx = 1:length(trials)
        
        % Done gathering spikes over trials
        
        if isempty(trialinfo)
            if norefOKflag
                warning('lfp_spikeAnalysis:nodata3', ...
                    'empty trialinfo table - implies no ref evts found');
                result = NaN;
                spikes = {};
                ntrigs = 0;
                return
            else
                error('lfp_spikeAnalysis:nodata2', ...
                    'empty trialinfo table - implies no ref evts found');
            end
        end
        ntrigs = size(trialinfo,1);
        % delete references to trials that were skipped because of no ref:
        spikes(norefidx,:) = [];
        if rasflag || placeflag || ~isempty(evts2avg)
            evtidx(norefidx,:) = [];
        end
        trialsincluded = trials(setdiff(1:length(trials), norefidx));
        
        if multiplotflag
            % may want to eventually account for more than 6 events...
            if fighandles(cluster) == 0
                if figflag
                    fignum = figbase + cluster;
                    if multichannelflag
                        fighandles(clustnums) = figure(fignum);
                    else
                        fighandles(cluster) = figure(fignum);
                    end
                else
                    if multichannelflag
                        fighandles(clustnums) = figure;
                    else
                        fighandles(cluster) = figure;
                    end
                end
            end
            % need this figtitle to be much shorter here b/c it labels each
            % subplot individually
            [pathstr, name] = fileparts(lfp_DataDir);
            figtitle = sprintf('%s\\%s\\%s align=%d', ...
                pathstr(max(1,end-6):end), name,...
                lfp_SpikeNames{cluster}, alignevents{i}(1)); %#ok<USENS>
            if length(alignevents{i}) > 1
                figtitle = [figtitle '...']; %#ok<AGROW>
            end
            figure(fighandles(cluster));
            if alignflag
                subplotdim1 = ceil(sqrt(length(alignevents)));
                subplotdim2 = floor(sqrt(length(alignevents)));
                if length(alignevents) > subplotdim1 * subplotdim2
                    subplotdim2 = subplotdim2 + 1;
                end
                hF = subplot(subplotdim1,subplotdim2,i);
                grid on;
            elseif multichannelflag
                hF = subplot(length(clustnums),1,clustidx);
                grid off;
            else
                grid on;
            end
        elseif meanfrflag
            % do nothing, just constructing mean firing rate between evts
            if nargout == 0 || returnfigflag
                warning('lfp_spikeAnalysis:nofig', ...
                    'type ''meanfr'' does not return a figure');
            end
        else
            if nargout == 0 || returnfigflag
                if figflag
                    fignum = figbase + cluster;
                    % close any previous incarnation of the figure:
                    try %#ok<TRYNC>
                        close(fignum);
                    end
                    hF = figure(fignum);
                else
                    if subplotsflag
                        hF = subplot(numrowscols(1), ...
                            numrowscols(2),plotnum);
                    else hF = figure;
                    end
                end
                if length(lfp_AlignmentRef) < 3
                    alignrefstr = mat2str(lfp_AlignmentRef);
                else
                    alignrefstr = sprintf('[%d...]', lfp_AlignmentRef(1));
                end
                figtitle = sprintf('%s\\%s align=%s win = %s', ...
                    lfp_DataDir, lfp_SpikeNames{cluster}, ...
                    alignrefstr, mat2str(window) );
                grid on;
            end
        end
        signifstr = '';
        
        % Do the maxskew test.
        % Note that it is logically impossible NOT to satisfy the condition
        % for returning NaN unless spikesfrac > numtrials.
        if percentage && nargout
            if multitrigflag
                error('lfp_spikeAnalysis:badopt3', ...
                    '''maxskew'' and ''multitrig'' cannot be combined' );
            end
            for trialindex = 1:size(spikes, 1)
                spikespertrial(trialindex) = length(spikes{trialindex}); %#ok<AGROW>
            end
            spikesfrac = round(sum(spikespertrial) * percentage/100);
            numtrials = round(trialfrac * length(spikes));
            spikespertrial = sort(spikespertrial);
            if sum(spikespertrial(end - numtrials + 1 : end)) > spikesfrac
                result = NaN;
                return
            end
        end

        if (rasflag || placeflag) && (nargout == 0 || returnfigflag)
            % Plot the raster or place field
            hold on;
            totalspikes = 0;
            if length(trialsincluded) ~= length(spikes)
                error('lfp_spikeAnalysis:oops', ...
                    'Dan screwed this one up!' );
            end
            for trialindex = 1:length(spikes)
                totalspikes = totalspikes + length(spikes{trialindex});
                if placeflag
                    lfp_plotplace(spikes{trialindex}, xchan, ychan);
                    hA = findobj(hF, 'Type', 'axes');
                end
                if rasflag
                    % 'multitrig' is explicitly disallowed for 'ras'
                    for myevtix = evtidx{trialindex,1}'
                        evtid = lfp_Events(myevtix,2);
                        if evtid == 0
                            warning('lfp_spikenalysis:eventID0', ...
                                'Event number %d has ID 0', myevtix );
                        elseif (evtid <= length(lfp_SelectedEventIDs)) ...
                                && lfp_SelectedEventIDs(evtid) ...
                                || (evtid > length(lfp_EventColors)) ...
                                && ~isempty(bigevts) ...
                                && ismember(evtid, bigevtIDs) %#ok<USENS>
                            % Either the event ID is selected, or it is a
                            % specified big event ID.
                            eventcolor = '';
                            eventname = ''; % required for detailstr
                            if evtid <= length(lfp_EventNames) %#ok<USENS>
                                eventname = lfp_EventNames{evtid};
                            end
                            if evtid <= length(lfp_EventColors)
                                eventcolor = lfp_EventColors{evtid};
                            else
                                for k = 1:size(bigevts,1)
                                    if ismember(evtid, bigevts{k,1}')
                                        eventcolor = bigevts{k,2};
                                        break
                                    end
                                end
                            end
                            if isempty(eventcolor)
                                eventcolor = lfp_EventDefaultColor;
                            end
                            hL = plot(...
                                lfp_Events(myevtix,1) - trialinfo(trialindex, 3), ...
                                trialindex, ...
                                'Color', eventcolor, ...
                                'Marker', '.' );
                            if isequal(class(eventcolor), 'char')
                                eventcolorstr = eventcolor;
                            else
                                eventcolorstr = mat2str(eventcolor);
                            end
                            detailstr = sprintf( ...
                                '\\nTimestamp=%.6f\\nLineColor="%s"\\nEventName="%s"\\nEventID=%.0f (0x%X)', ...
                                lfp_Events(myevtix,1), ...
                                eventcolorstr, eventname, evtid, evtid );
                            set(hL, ...
                                'ButtonDownFcn', ...
                                ['fprintf(1,''' detailstr '\n'')'] );
                        end
                    end 
                    dg_plottick(spikes{trialindex}, trialindex, 1, 'k');
                end %if rasflag
            end %for trialindex = 1:length(spikes)
            if rasflag % #2
                hA = findobj(hF, 'Type', 'axes');
                set(hA,'YDir','reverse');
                set(get(hA,'YLabel'), 'String', 'Trial Number');
                set(get(hA,'XLabel'), 'String', 'Time, seconds');
                yticklabels = get(hA, 'YTickLabel');
                % All this yticklabel stuff gets very unhappy when using
                % subplots; this is just a band-aid on the festering sore:
                if iscell(yticklabels)
                    yticklabels = char(yticklabels);
                end
                % end band-aid
                newyticklabels = {};
                for ytidx = 1:size(yticklabels,1)
                    ytn = str2double(yticklabels(ytidx,:));
                    if (fix(ytn) == ytn) && (ytn ~= 0) && (ytn <= length(trialsincluded))
                        newyticklabels{ytidx} = num2str(trialsincluded(ytn)); %#ok<AGROW>
                    else
                        newyticklabels{ytidx} = ''; %#ok<AGROW>
                    end
                end
                set(hA, 'YTickLabel', newyticklabels);
                maxtimerange = [min(trialinfo(:,1)) max(trialinfo(:,2))];
                xlim(hA, maxtimerange);
            end %if rasflag % #2
        end

        if hisflag || cdfflag
            % Calculate and plot the histogram
            if any(trialinfo(:,3) == 0)
                error('lfp_spikeAnalysis:missingref', ...
                    'Some selected trials have no reference event.' );
            end
            starttime = max(trialinfo(:,1));
            endtime = min(trialinfo(:,2));
            dur = endtime - starttime;
            if newbinwidthflag || newbinsflag
                starttime = (round(1000000*starttime)/1000000);
                dur = (round(1000000*dur)/1000000);
                if rawbinwidth == 0
                    % nbins was specified; calculate max integral ms
                    % binwidth that will accomodate nbins.
                    binwidth = fix(dur*1000/nbins);
                    if binwidth == 0
                        close(hF);
                        error('lfp_spikeAnalysis:binwidth', ...
                            'The time window of %.6f s is too small to accommodate %d bins', ...
                            dur, nbins);
                    end
                else
                    binwidth = rawbinwidth;
                    nbins = fix(dur*1000/binwidth);
                end
                bins = starttime + ((1:nbins) - 1/2)*binwidth/1000;
                for cellidx = 1:numel(spikes)
                    spikes{cellidx} = ...
                        round(spikes{cellidx} * 1000000) / 1000000 ;
                end
                binwidth = binwidth * 1e-3; % finally convert to s from ms
            else
                binwidth = rawbinwidth;
                % calculate bin centers based on actual start and end time and
                % specified binwidth or nbins; if time 0 is included, it should
                % be a bin center, otherwise the first bin center should be at
                % starttime + 0.5*binwidth.
                if (starttime < 0) && (endtime > 0)
                    % time 0 is included
                    if rawbinwidth == 0
                        % nbins was specified; start with approximate value of
                        % binwidth to locate the bin that contains time zero,
                        % and then make binwidth small enough so that the first
                        % bin and last bin are both entirely contained within
                        % starttime and endtime.
                        binw1 = dur/(nbins);
                        nbinsbefore = floor(-starttime/binw1);
                        nbinsafter = nbins - nbinsbefore - 1;
                        binwidth = min(-starttime/(nbinsbefore + 0.5), ...
                            endtime/(nbinsafter + 0.5) );
                    else
                        % binwidth was specified
                        binwidth = rawbinwidth;
                        nbinsbefore = floor((-starttime-0.5*binwidth)/binwidth);
                        nbinsafter = floor((endtime-0.5*binwidth)/binwidth);
                        nbins = nbinsbefore + nbinsafter + 1;
                    end
                    bins = (-nbinsbefore : nbinsafter) * binwidth;
                else
                    % time 0 is not included
                    if binwidth == 0
                        % nbins was specified
                        binwidth = dur/(nbins);
                    else
                        % binwidth was specified
                        nbins = floor(dur/binwidth);
                    end
                    bins = starttime + (0.5:nbins)*binwidth;
                end

                if offsetflag
                    bins = bins - binwidth/2;
                end
            end
            
            % trim off spikes that go beyond the range of bins, and count the
            % remaining spikes.
            totalspikes = 0;
            for cellidx = 1:numel(spikes)
                spikes{cellidx} = spikes{cellidx}( ...
                     spikes{cellidx} >= bins(1) - 0.5 * binwidth ...
                    & spikes{cellidx} <= bins(end) + 0.5 * binwidth );
                totalspikes = totalspikes + length(spikes{cellidx});
            end
            mergedspikes = zeros(totalspikes,1);
            nmerged = 0;
            for cellidx = 1:numel(spikes)
                mergedspikes( ...
                    nmerged + 1 : nmerged + numel(spikes{cellidx}) ...
                ) = spikes{cellidx} ;
                nmerged = nmerged + numel(spikes{cellidx});
            end
            if hzflag
                % Number of rows in trialinfo reflects the actual number of
                % triggers whether multitrig or not:
                scale = 1/(binwidth * ntrigs);
            else
                scale = 1;
            end
            if nargout > 0 && ~returnfigflag
                if cdfflag
                    result = reshape(cdfplot(mergedspikes), [], 1);
                else
                    result = reshape( ...
                        dg_hist(mergedspikes, bins, scale), [], 1 );
                end
            else
                if cdfflag
                    if ~isempty(mergedspikes)
                        cdfplot(mergedspikes);
                    end
                else
                    dg_hist(mergedspikes, bins, scale);
                end
                hA = findobj(hF, 'Type', 'axes');
                if trigfuncflag || isempty(lfp_AlignmentRef)
                    lfp_plotEvtMarkers(hA, [0 NaN]);
                else
                    lfp_plotEvtMarkers(hA, [0 lfp_AlignmentRef(1)]);
                end
                xlabel('Time, seconds');
                if (subplotsflag || multiplotflag) && multitrigflag
                    if hzflag
                        ylabel(sprintf('Hz, %d/%d', ...
                            ntrigs, length(trialsincluded)));
                    else
                        ylabel(sprintf('Spikes, %d/%d', ...
                            ntrigs, length(trialsincluded)));
                    end
                else
                    if hzflag
                        ylabel(sprintf('Hz, %d trials', ...
                            length(trialsincluded)));
                    else
                        ylabel(sprintf('Spikes, %d trials', ...
                            length(trialsincluded)));
                    end
                end
            end
        end %if hisflag || cdfflag
        
        % For phase histograms, the method is to find positive-going zero
        % crossings in the reference waveform in wavefilenum, then for each
        % interval between successive zero crossings, add a count for each
        % spike to the appropriate phase bin.
        if phaseflag
            % Note that if offsetflag is not set, then we get an extra
            % bincenter because "0" and "1" bins are only a half-bin wide.
            bincenters = ( offsetflag*0.5 : nbins )/nbins;
            histos = zeros(length(trials), length(bincenters));
            for trialindex = 1:length(trials)
                timerange = trialinfo(trialindex, 1:2) + trialinfo(trialindex, 3);
                samplerange(1) = lfp_time2index(timerange(1));
                samplerange(2) = lfp_time2index(timerange(2));
                
                % find zero crossing timestamps (xingTS):
                pos_xings = lfp_findThresholds(0, 0, ...
                    wavefilenum, samplerange);
                pos_xings = pos_xings + samplerange(1);
                xingTS = zeros(1, length(pos_xings));
                for xingindex = 1:length(pos_xings)
                    xingTS(xingindex) = lfp_index2time(pos_xings(xingindex));
                end
                % make xingTS relative to reftime, like spikes
                xingTS = xingTS - trialinfo(trialindex, 3);
                
                % Construct phase histo for this trial; norm_spikes is
                % normalized such that one cycle = 1.0
                for tsindex = 2:length(xingTS)
                    ts1 = xingTS(tsindex - 1);
                    ts2 = xingTS(tsindex);
                    onecyclespikes = spikes{trialindex}( ...
                        spikes{trialindex} >= ts1 & spikes{trialindex} < ts2 ...
                        );
                    norm_spikes = (onecyclespikes - ts1) / (ts2 - ts1);
                    histos(trialindex, :) = ...
                        histos(trialindex, :) + hist(norm_spikes, bincenters);
                end
            end
            sumhisto = sum(histos, 1);
            totalspikes = sum(sumhisto);
            signif = dg_chi2test2(sumhisto, true);
            signifstr = sprintf(' p<%.4f', signif);
            if ~offsetflag
                % The bin centered on zero must be recalculated by adding the
                % counts in the raw "0" bin, which is only half a binwidth
                % wide, to the counts in the raw "1" bin, which is also only
                % half a binwidth wide.  Then we throw away the raw "1" bin.
                sumhisto(1) = sumhisto(1) + sumhisto(end);
                sumhisto(end) = [];
                bincenters(end) = [];
            end
            meanspikes = totalspikes/length(sumhisto);
            sd = std(sumhisto);
            angles = 2*pi*[bincenters bincenters(1)];
            pluserr = meanspikes+2*sd;
            minuserr = meanspikes-2*sd;
            if nargout > 0
                result = reshape(sumhisto, [], 1);
            else
                hF = figure(hF);
                % include pluserr to set scale
                polar([angles' angles'], ...
                    [repmat(pluserr, size(angles))' ...
                        [sumhisto sumhisto(1)]' ...
                    ], 'k');
                hold on;
                % draw pluserr and minuserr in correct colors
                polar(angles', repmat(pluserr, size(angles))', 'r');
                if minuserr < 0
                    warning('lfp_spikeAnalysis:minuserrLT0', ...
                        'Omitting mean - 2 SD line due to value less than zero.');
                else
                    polar(angles', repmat(minuserr, size(angles))', 'b');
                end
                hA = findobj(hF, 'Type', 'axes');
                hline = findobj(hA,'Type','line');
                set(hline,'LineWidth',1.5);
                set(hline(end - 1), 'Marker', 'o');
                figtitle = [ figtitle ' ' lfp_FileNames{wavefilenum} ]; %#ok<AGROW>
            end
        end
        
        if nargout == 0 || returnfigflag
            if ~multiplotflag
                if SeqTypeTitleflag
                %this adds spatial and temporal parameters of the sequence
%                     figtitle = sprintf(...
%                         '%s N=%d%s\n%s %d%d%d %d%d%d', ...
%                         figtitle, totalspikes, signifstr, ...
%                         lfp_getTrialsLabel(trials, trialstyle), ...
%                         lfp_TrialParams{trials(1)}(8:10), ...
%                         lfp_TrialParams{trials(1)}(12:14));
                %this adds only temporal parameters of the sequence
                    figtitle = sprintf(...
                        '%s N=%d%s\n%s XXX %d%d%d', ...
                        figtitle, totalspikes, signifstr, ...
                        lfp_getTrialsLabel(trials, trialstyle), ...
                        lfp_TrialParams{trials(1)}(12:14)); %#ok<USENS>
                else
                    figtitle = sprintf(...
                        '%s N=%d%s\n%s', ...
                        figtitle, totalspikes, signifstr, ...
                        lfp_getTrialsLabel(trials, trialstyle));
                end
                if hisflag
                    if multitrigflag
                        trigcount = ntrigs;
                    else
                        trigcount = 0;
                    end
                    if newbinwidthflag || newbinsflag
                        figtitle = sprintf(...
                            '%s newbins=%d binwidth=%d multitrig=%d', ...
                            figtitle, nbins, 1e+3*binwidth, trigcount );
                    else
                        figtitle = sprintf(...
                            '%s nbins=%d binwidth=%.2f multitrig=%d', ...
                            figtitle, nbins, 1e+3*binwidth, trigcount );
                    end
                    if trigfuncflag
                        figtitle = sprintf('%s %s', ...
                            figtitle, func2str(trigfuncH));
                    end
                elseif phaseflag
                    figtitle = sprintf(...
                    '%s nbins=%d', ...
                    figtitle, nbins );
                end
            end
            hT = get(hA, 'Title');
            if titleflag
                figtitle = titlestr;
            end
            set(hT, 'Interpreter', 'none', 'String', figtitle);
        end
        
        if ~isempty(evts2avg) && ~(placeflag || phaseflag)
            % Draw aggregate events markers.
            % <evtmatrix> contains events in lfp_Events format, but with
            % time relative to the reference event for each trigger.
            hA = findobj(hF, 'Type', 'axes');
            lfp_plotEvtMarkers(hA, evtmatrix, 'stats', 'bigevts', bigevts);
        end %if ~isempty(evts2avg) && ~(placeflag || phaseflag)
        
        if meanfrflag
            % This requires no figure and only sets result to be the mean
            % firing rate (hz) over all time
            totalspikes = length(cat(1,spikes{:}));
            totaltime = sum(trialinfo(:,2)-trialinfo(:,1));
            
            pertrialtimes = trialinfo(:,2)-trialinfo(:,1);
            pertrialspikes = zeros(size(trialinfo,1),1);
            pertrialmean = zeros(size(trialinfo,1),1);
            for k = 1:size(trialinfo,1)
                pertrialspikes(k) = length(spikes{k});
                pertrialmean(k) = pertrialspikes(k)/pertrialtimes(k);
            end            
            result = {totalspikes/totaltime pertrialmean};
        end
        
    end %for clustidx = 1:length(clustnums)
end %for i = 1:length(alignevents)

if alignflag
    lfp_AlignmentRef = oldAlignmentRef; %#ok<NASGU>
end
