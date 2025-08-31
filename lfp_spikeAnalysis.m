function [result, spikes, ntrigs, hF] = lfp_spikeAnalysis( ...
    type, trials, clustnums, window, varargin)
%LFP_SPIKEANALYSIS provides histograms (PETH and phase) and rasters.
%lfp_spikeAnalysis(type, trials, clustnums, window)
%lfp_spikeAnalysis(type, trials, clustnums)
%lfp_spikeAnalysis(type, trials)
%lfp_spikeAnalysis(type)
%lfp_spikeAnalysis(..., 'bigdata')
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
%lfp_spikeAnalysis(..., 'returntimepts')
%lfp_spikeAnalysis(..., 'SeqTypeTitle')
%lfp_spikeAnalysis(..., 'subplots', numrowscols, plotnum)
%lfp_spikeAnalysis(..., 'title', titlestr)
%lfp_ras(...)
%lfp_ras(..., 'ticksize', width)
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
%  <spikes> is a cell array of column vectors of spike times relative to
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
%   'meanfr' - can be called with 'evtbounds' option, or with the standard
%       alignment and window combination. result{1} is the total
%       spikes/total time and result{2} is a vector containing the mean
%       firing rate per trial.  Do not specify more than one cluster at a
%       time together with this type of analysis.  The 'hz' flag does not
%       change the result of a 'meanfr' analysis.
%   'phase' - invoked by lfp_phase(..., wavefilenum); <result> is a column
%       vector. For phase histograms, the method is to find positive-going
%       zero crossings in the reference waveform in wavefilenum, then for
%       each interval between successive zero crossings, add a count for
%       each spike to the appropriate phase bin.
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
%       greater than 0.
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
%       (Note that this is exactly the opposite to the behavior of 'offset'
%       in lfp_JPSTH.  Sorry about that!  But it's too late now...)
%   If <type> is 'phase', then the 'nbins' option works as for 'his'.
%lfp_spikeAnalysis(..., 'aligns', alignmentlist)
%  Creates up to six subplots for up to six different alignment events in
%  one figure.  Each of the event IDs in the cell array <alignmentlist> is
%  used in turn as the alignment event; this is done by actually assigning
%  the value to lfp_AlignmentRef, so any previous value of lfp_AlignmentRef
%  is lost.
%lfp_spikeAnalysis(..., 'bigdata')
%  Skips the error check for too many spikes in raster plot.  Note that
%  if the number of spikes is large enough, Matlab could freeze up your
%  whole machine.
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
%  lfp_XLimAll and <window> are both empty.  {This code looks like it
%  cannot work if more than one cluster is specified!  Also note that
%  'multitrig' is allowed ONLY when <type> = 'his'. DG 26-Aug-2013}
%lfp_spikeAnalysis(..., 'ticksize', width)
%  This only affects 'ras' plots, and sets the 'LineWidth' property of the
%  spike tickmarks to <width>.
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
%  Forces creation of a figure when return values are used; if 'returnfig'
%  is not used, then hF = [].
%lfp_spikeAnalysis(..., 'returntimepts')
%  When <type> is 'his', this causes <result> to be a two-column array
%  where the first column contains bin center times.
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

% Regression test: test_lfp_spikeAnalysis.

%$Rev: 421 $
%$Date: 2023-03-31 20:31:31 -0400 (Fri, 31 Mar 2023) $
%$Author: dgibson $

global lfp_XLimAll lfp_Spikes lfp_TrialIndex lfp_FileNames lfp_TrialStyle
global lfp_AlignmentRef lfp_Events lfp_DataDir lfp_SpikeNames
global lfp_SelectedEventIDs lfp_EventColors lfp_EventNames
global lfp_EventDefaultColor lfp_EventShapes lfp_EventDefaultShape
global lfp_TrialParams

result = [];
spikes = [];
ntrigs = [];
hF = [];

if nargin <4 || isempty(window)
    window = lfp_XLimAll;
elseif ~(isa(window, 'double') ...
        && isequal(size(window), [1 2]) )
    error('lfp_spikeAnalysis:badwindow', ...
        '<window> must be 1x2 number array.' );
end
if nargin < 3 || isempty(clustnums)
    clustnums = 1:length(lfp_Spikes);
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
gatherTrialSpikesOpts = {}; % subsequent code guarantees column vector
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
        gatherTrialSpikesOpts{end+1,1} = 'place'; %#ok<*NASGU>
    case 'ras'
        rasflag = true;
    case 'meanfr'
        meanfrflag = true;
        gatherTrialSpikesOpts{end+1,1} = 'meanfr';
    otherwise
        error('lfp_spikeAnalysis:badMode', ...
            '<type> must be one of: ras, his.' );
end

if phaseflag
    if ~( (~isempty(varargin)) ...
            && isa(varargin{1}, 'double') ...
            && (fix(varargin{1}) == varargin{1}) )
        error('lfp_spikeAnalysis:needfilenum', ...
            'You must specify a wave filenum after <window> to use ''phase''' );
    end
    wavefilenum = varargin{1};
    varargin(1) = [];
end

if placeflag
    xchan = find(ismember(lfp_FileNames, 'trk X'));
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
bigdataflag = false;
bigevts = {};
evtavg2flag = false;
evtbounds = {};
evts2avg = [];
% The <figflag>, <returnfigflag>, <plotflag> rat's nest: they are both
% false by default. figflag becomes true if a figure number is specified
% (horrible feature that nobody uses and causes undue complexities).
% returnfigflag is totally unrelated, and controls figure window creation.
% The condition for creating a figure is (nargout == 0 || returnfigflag).
% All calls to 'figure' that create a new figure assign the figure handle
% to either hF or the appropriate element of fighandles.  However, the
% calls to plotting routines ('plot', 'polar', 'bar', 'dg_plottick',
% 'lfp_plotEvtMarkers', etc) do not reference any graphics handles (this
% should also be fixed), so that means they can also create figures. Since
% we don't know a priori whether to look at hF or fighandles to see if a
% figure has been created, all calls to plotting routines need to be
% withing the scope of an "if nargout == 0 || returnfigflag".  I have
% therefore provided the syntactic sugar of plotflag = (nargout == 0 ||
% returnfigflag) to make it easy to more-or-less immediately surround every
% plotting call in an "if plotflag".  Clearly, major refactoring is in
% order someday...
figflag = false;
returnfigflag = false;
norefOKflag = false;
offsetflag = false;
alignflag = false;
multichannelflag = false;
multitrigflag = false;
newbinsflag = false;
newbinwidthflag = false;
SeqTypeTitleflag = false;
nbins = 100;
rawbinwidth = 0;
alignevents = {}; %TF, DG
hzflag = false; %TF 2/21/05
returntimeptsflag = false;
tickwidth = 0.5;
trigfuncflag = false;
titleflag = false;
while argnum <= length(varargin)
    if isa(varargin{argnum}, 'double')
        if figflag
            error('lfp_spikeAnalysis:multifigbase', ...
                'You have specified figbase more than once.');
        else
            figflag = true;
            figbase = varargin{argnum};
        end
    else
        switch varargin{argnum}
            case 'aligns' %TF added, DG modified
                argnum = argnum +1;
                alignflag = true;
                alignevents = varargin{argnum};
                oldAlignmentRef = lfp_AlignmentRef;
            case 'bigdata'
                bigdataflag = true;
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
                        || (~isa(varargin{argnum}, 'double'))
                    error('lfp_spikeAnalysis:badbinwidth', ...
                        'You must specify a numeric value for binwidth' );
                end
                rawbinwidth = 1.0e-3 * varargin{argnum}; % convert ms to s
            case 'evtavg'
                argnum = argnum + 1;
                evts2avg = varargin{argnum};
                gatherTrialSpikesOpts{end+1,1} = 'evtavg'; %#ok<*AGROW>
                gatherTrialSpikesOpts{end+1,1} = evts2avg;
            case 'evtavg2'
                argnum = argnum + 1;
                evts2avg = varargin{argnum};
                evtavg2flag = true;
                gatherTrialSpikesOpts{end+1,1} = 'evtavg2';
                gatherTrialSpikesOpts{end+1,1} = evts2avg;
            case 'evtbounds'
                argnum = argnum + 1;
                evtbounds = varargin{argnum};
                gatherTrialSpikesOpts{end+1,1} = 'evtbounds';
                gatherTrialSpikesOpts{end+1,1} = evtbounds;
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
                gatherTrialSpikesOpts{end+1,1} = 'multitrig';
            case 'nbins'
                argnum = argnum + 1;
                if (argnum > length(varargin)) ...
                        || (~isa(varargin{argnum}, 'double'))
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
                gatherTrialSpikesOpts{end+1,1} = 'norefOK';
            case 'offset'
                offsetflag = true;
            case 'returnfig'
                returnfigflag = true;
            case 'returntimepts'
                returntimeptsflag = true;
            case 'SeqTypeTitle'
                SeqTypeTitleflag = true;
            case 'session'
                argnum = argnum + 1;
                session = varargin{argnum};
            case 'showtrialnums'
                trialstyle = 'trialnums';
            case 'subplots'
                subplotsflag = true;
                argnum = argnum + 1;
                numrowscols = varargin{argnum};
                argnum = argnum + 1;
                plotnum = varargin{argnum};
            case 'ticksize'
                argnum = argnum + 1;
                tickwidth = varargin{argnum};
            case 'title'
                titleflag = true;
                argnum = argnum + 1;
                titlestr = varargin{argnum};
            case 'trigfunc'
                trigfuncflag = true;
                argnum = argnum +1;
                trigfuncH = varargin{argnum};
                trigfuncArgs = varargin(argnum+1:end);
                argnum = length(varargin);
                gatherTrialSpikesOpts = [
                    gatherTrialSpikesOpts
                    {'trigfunc'}
                    {trigfuncH}
                    reshape(trigfuncArgs, [], 1)
                    ];
            otherwise
                error('lfp_spikeAnalysis:badoption', ...
                    ['The option "' varargin{argnum} '" is not recognized.'] );
        end
    end
    argnum = argnum + 1;
end

if meanfrflag
    if isempty(evtbounds) && isempty(window)
        error('lfp_spikeAnalysis:betweenevtbounds', ...
            '''meanfr'' requires either ''evtbounds'' or a non-empty <window>');
    end
    if  ~isempty(evtbounds) && ~isempty(window)
        error('lfp_spikeAnalysis:betweenevtwindow', ...
            '''meanfr'' with ''evtbounds'' requires empty <window>');
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

if nargout ~= 0 && ~returnfigflag && ~isempty(evts2avg)
    warning('lfp_spikeAnalysis:nofig', ...
        'Ignoring ''evts2avg'' option because no figure is being created');
    evts2avg = [];
end

if ~isa(nbins, 'double')
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
if ~(isa(trials, 'double')) ...
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

plotflag = (nargout == 0 || returnfigflag);

% run the whole procedure for each alignment ref
for i = 1:length(alignevents)
    if length(alignevents) > 1 || multichannelflag
        multiplotflag = true;
    else
        multiplotflag = false;
    end
    lfp_AlignmentRef = alignevents{i};% DG modified
    
    
    %%% Large block of code excised and replace with: %%%%%%%%%%%%%%%%%
    if rasflag || placeflag || ~isempty(evts2avg)
        if trigfuncflag
            for argnum = 1:length(gatherTrialSpikesOpts)
                if isequal(gatherTrialSpikesOpts{argnum}, 'trigfunc')
                    break
                end
            end
            gatherTrialSpikesOpts = [
                gatherTrialSpikesOpts(1:argnum-1)
                {'evts'}
                gatherTrialSpikesOpts(argnum:end)
                ];
        else
            gatherTrialSpikesOpts{end+1,1} = 'evts';
        end
    end
    [spikes, trialinfo, evtidx, evtmatrix, noreftrialidx] = ...
        lfp_getSpikes( trials, clustnums, window, ...
        gatherTrialSpikesOpts{:} );
    %%% End replacement for excised block of code. %%%%%%%%%%%%%%%%%%%%
    
    for clustidx = 1:length(clustnums)
        clusternum = clustnums(clustidx);
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
        spikes(noreftrialidx, :, :) = [];
        if rasflag || placeflag || ~isempty(evts2avg)
            evtidx(noreftrialidx,:) = [];
        end
        trialsincluded = trials(setdiff(1:length(trials), noreftrialidx));
        
        if multiplotflag
            % may want to eventually account for more than 6 events...
            if fighandles(clusternum) == 0
                if figflag
                    fignum = figbase + clusternum;
                    if multichannelflag
                        fighandles(clustnums) = figure(fignum);
                    else
                        fighandles(clusternum) = figure(fignum);
                    end
                else
                    if multichannelflag
                        fighandles(clustnums) = figure;
                    else
                        fighandles(clusternum) = figure;
                    end
                end
            end
            % need this figtitle to be much shorter here b/c it labels each
            % subplot individually
            [pathstr, name] = fileparts(lfp_DataDir);
            figtitle = sprintf('%s\\%s\\%s align=%d', ...
                pathstr(max(1,end-6):end), name,...
                lfp_SpikeNames{clusternum}, alignevents{i}(1));
            if length(alignevents{i}) > 1
                figtitle = [figtitle '...']; 
            end
            figure(fighandles(clusternum));
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
            if plotflag
                warning('lfp_spikeAnalysis:nofig', ...
                    'type ''meanfr'' does not return a figure');
            end
        else
            if plotflag
                if figflag
                    fignum = figbase + clusternum;
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
                    lfp_DataDir, lfp_SpikeNames{clusternum}, ...
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
                spikespertrial(trialindex) = length(spikes{trialindex}); 
            end
            spikesfrac = round(sum(spikespertrial) * percentage/100);
            numtrials = round(trialfrac * size(spikes,1));
            spikespertrial = sort(spikespertrial);
            if sum(spikespertrial(end - numtrials + 1 : end)) > spikesfrac
                result = NaN;
                return
            end
        end
        
        if (rasflag || placeflag) && plotflag
            totalspikes = numel(cell2mat(reshape( ...
                spikes(:,1,clustidx), [], 1 )));
            if rasflag && ~bigdataflag && totalspikes > 1e4
                error('lfp_spikeAnalysis:bigdata', ...
                    'Too many spikes (%d) for raster plot; please reduce time window or number of trials.', ...
                    totalspikes);
            end
            % Plot the raster or place field
            hold on;
            if length(trialsincluded) ~= size(spikes,1)
                error('lfp_spikeAnalysis:oops', ...
                    'Dan screwed this one up!' );
            end
            for trialindex = 1:size(spikes,1)
                totalspikes = totalspikes + length(spikes{trialindex,1,clustidx});
                if placeflag
                    lfp_plotplace(spikes{trialindex,1,clustidx}, xchan, ychan);
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
                                && ismember(evtid, bigevtIDs)
                            % Either the event ID is selected, or it is a
                            % specified big event ID.
                            eventcolor = '';
                            eventname = ''; % required for detailstr
                            if evtid <= length(lfp_EventNames)
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
                            if evtid <= length(lfp_EventShapes) ...
                                    && ~isempty(lfp_EventShapes{evtid})
                                eventshape = lfp_EventShapes{evtid};
                            else
                                eventshape = lfp_EventDefaultShape;
                            end
                            if plotflag
                                hL = plot(...
                                    lfp_Events(myevtix,1) - trialinfo(trialindex, 3), ...
                                    trialindex, ...
                                    'Color', eventcolor, ...
                                    'MarkerSize', 12, ...
                                    'Marker', eventshape );
                            end
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
                    if plotflag
                        dg_plottick(spikes{trialindex,1,clustidx}, ...
                            trialindex, 1, 'k', tickwidth);
                    end
                end %if rasflag
            end %for trialindex = 1:size(spikes,1)
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
                        newyticklabels{ytidx} = num2str(trialsincluded(ytn)); 
                    else
                        newyticklabels{ytidx} = ''; 
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
                if offsetflag
                    [binedges, bins] = dg_makeTimeBins( ...
                        [starttime, endtime], binwidth);
                else
                    [binedges, bins] = dg_makeTimeBins( ...
                        [starttime - binwidth/2000, endtime], binwidth); %#ok<*ASGLU>
                end
                % round to integral microseconds
                for trialidx = 1:size(spikes,1)
                    for trigidx = 1:size(spikes,2)
                        spikes{trialidx, trigidx, clustidx} = ...
                            round( spikes{trialidx, trigidx, clustidx} ...
                            * 1000000 ) / 1000000 ;
                    end
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
            for trialidx = 1:size(spikes,1)
                for trigidx = 1:size(spikes,2)
                    spikes{trialidx, trigidx, clustidx} = ...
                        spikes{trialidx, trigidx, clustidx}( ...
                        spikes{trialidx, trigidx, clustidx} >= bins(1) - 0.5 * binwidth ...
                        & spikes{trialidx, trigidx, clustidx} <= bins(end) + 0.5 * binwidth );
                    totalspikes = totalspikes + ...
                        length(spikes{trialidx, trigidx, clustidx});
                end
            end
            mergedspikes = zeros(totalspikes,1);
            nmerged = 0;
            for trialidx = 1:size(spikes,1)
                for trigidx = 1:size(spikes,2)
                    mergedspikes( ...
                        nmerged + 1 : nmerged + ...
                        numel(spikes{trialidx, trigidx, clustidx}) ...
                        ) = spikes{trialidx, trigidx, clustidx} ;
                    nmerged = nmerged + numel(spikes{ ...
                        trialidx, trigidx, clustidx});
                end
            end
            if hzflag
                % Number of rows in trialinfo reflects the actual number of
                % triggers whether multitrig or not:
                scale = 1/(binwidth * ntrigs);
            else
                scale = 1;
            end
            if cdfflag
                if isempty(mergedspikes)
                    result = 0;
                else
                    result = reshape(cdfplot(mergedspikes), [], 1);
                end
            else
                [n, x] = dg_hist(mergedspikes, bins, scale);
                if plotflag
                    bar(x, n, 'hist');
                end
                if returntimeptsflag
                    result = [reshape(x, [], 1) reshape(n, [], 1)];
                else
                    result = reshape(n, [], 1);
                end
            end
            if plotflag
                hA = findobj(hF, 'Type', 'axes');
                if trigfuncflag || isempty(lfp_AlignmentRef)
                    lfp_plotEvtMarkers(hA, [0 NaN]);
                else
                    lfp_plotEvtMarkers(hA, [0 lfp_AlignmentRef(1)], ...
                        'bigevts', bigevts);
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
        
        if phaseflag
            % Note that if offsetflag is not set, then we get an extra
            % bincenter because "0" and "1" bins are only a half-bin wide.
            % This issue is handled later.
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
                    onecyclespikes = spikes{trialindex, 1, clustidx}( ...
                        spikes{trialindex, 1, clustidx} >= ts1 & spikes{trialindex, 1, clustidx} < ts2 ...
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
            end
            if plotflag
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
                figtitle = [ figtitle ' ' lfp_FileNames{wavefilenum} ]; 
            end
        end
        
        if plotflag
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
                        lfp_TrialParams{trials(1)}(12:14));
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
        
        if plotflag && (~isempty(evts2avg) && ~(placeflag || phaseflag))
            % Draw aggregate events markers.
            % <evtmatrix> contains events in lfp_Events format, but with
            % time relative to the reference event for each trigger.
            hA = findobj(hF, 'Type', 'axes');
            lfp_plotEvtMarkers(hA, evtmatrix, 'stats', 'bigevts', bigevts);
        end %if ~isempty(evts2avg) && ~(placeflag || phaseflag)
        
        if meanfrflag
            % It matters not whether <spikes> was selected per 'evtbounds'
            % option or per align/window method.  The computation remains
            % the same:
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
    lfp_AlignmentRef = oldAlignmentRef;
end
