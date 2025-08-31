function [hF, Spectra, f, phi, ntrigs, Spectra2plot] = lfp_mtspectrum(...
    trials, filenums, window, varargin)
%LFP_MTSPECTRUM shows power spectra and coherence using chronux functions.
%lfp_mtspectrum(trials, filenums, window)
%lfp_mtspectrum(..., 'arrows', n)
%lfp_mtspectrum(..., 'avg')
%lfp_mtspectrum(..., 'coh')
%lfp_mtspectrum(..., 'EP')
%lfp_mtspectrum(..., 'err', [n p])
%lfp_mtspectrum(..., figbase)
%lfp_mtspectrum(..., 'k', numtapers)
%lfp_mtspectrum(..., 'lines', n)
%lfp_mtspectrum(..., 'log')
%lfp_mtspectrum(..., 'multiwin', winwidth)
%lfp_mtspectrum(..., 'nodisplay')
%lfp_mtspectrum(..., 'noplot')
%lfp_mtspectrum(..., 'norm')
%lfp_mtspectrum(..., 'nw', nw)
%lfp_mtspectrum(..., 'ovr')
%lfp_mtspectrum(..., 'p', plevel)
%lfp_mtspectrum(..., 'pad', N)
%lfp_mtspectrum(..., 'phi')
%lfp_mtspectrum(..., 'print')
%lfp_mtspectrum(..., 'rmBL', BL)
%lfp_mtspectrum(..., 'rmdc')
%lfp_mtspectrum(..., 'rmEP')
%lfp_mtspectrum(..., 'rotate')
%lfp_mtspectrum(..., 'session', sessionname)
%lfp_mtspectrum(..., 'showtrialnums')
%lfp_mtspectrum(..., 'thresh')
%lfp_mtspectrum(..., 'unwrap')
%lfp_mtspectrum(..., 'xscale', scale)
%lfp_mtspectrum(..., 'yscale', scale)
%lfp_mtspectrum(trials, filenums)
%lfp_mtspectrum(trials)
%lfp_mtspectrum

%lfp_mtspectrum(trials, filenums, window)
%  <trials> is as in lfp_disp
%  Plots the FFT computed on time window <window> (which is a 2-element
%  array specifying [start-time end-time]) from the file numbers in
%  <filenums>.  If <filenums> is the empty array or not given, then all
%  files are plotted.  If <window> is empty or not given, it defaults to
%  lfp_XLimAll, and if that is also empty is defaults to [0 1].  For each
%  trial, <window> is clipped to the times of the lfp_NominalTrialStart and
%  lfp_NominalTrialEnd events if lfp_XLimAll is empty, or to the times of
%  the beginning and end of the surrounding recorded segment otherwise.
%  Note that lfp_XLimAll still controls how the window gets clipped even
%  when <window> is not empty. Time is shown relative to the first
%  lfp_AlignmentRef event between the trial's lfp_NominalTrialStart and
%  lfp_NominalTrialEnd events.  It is an error if there is no
%  lfp_AlignmentRef event in one of the trials.
%lfp_mtspectrum(..., 'arrows', n)
%  Works only with 'coh'.  Draws little phase arrows on the curve in the
%  regions where it is above the p<.01 threshold.  <n> specifies how many
%  frequency points at which to test threshold and potentially plot arrows.
%  When the arrows are pointing clockwise of zero (zero being horizontal to
%  the right, and clockwise being negative phase), that means that the
%  second channel is lagging the first.
%lfp_mtspectrum(..., 'avg')
%  Instead of plotting each trial in its own figure window, the data for
%  the trials are averaged after being FFT'd, and then plotted.
%lfp_mtspectrum(..., 'coh')
%  Computes the coherence spectrum between two channels. <filenums> must
%  have an even number of entries after being masked by lfp_SelectedFiles;
%  a separate plot will be shown for each pair. Vertical scale is [0 1]
%  without error bars or [0 1.5] with error bars.
%lfp_mtspectrum(..., 'EP')
%  Averages the waveforms first, then computes spectrum of the average.
%lfp_mtspectrum(..., 'err', [n p])
%  Displays chronux error curves for significance level <p>; <n>=1 for
%  theoretical calculation, <n>=2 for jackknife (takes longer).  When
%  combined with 'coh', <n>=1 just adds a horizontal line at the coherence
%  threshold for statistically significant difference from 0.
%lfp_mtspectrum(..., 'evtbounds', {startIDs endIDs})
%  Sets the range of time to search for alignment events to something
%  narrower than the entire trial, specified in terms of other events.
%  <startIDs> is a list of alternative event IDs to use as the start of the
%  time range, and <endIDs> is a list of alternative event IDs to use as
%  the end of the time range.  The start event is the earliest event in the
%  trial that is a member of <startIDs>; the end event is the earliest
%  event in the trial AFTER the start event that is a member of <endIDs>.
%lfp_mtspectrum(..., 'k', numtapers)
%  Forces number of multitapers to <numtapers>; default is 1.
%lfp_mtspectrum(..., 'lines', n)
%  Same as 'arrows', except it just draws little line segments (good for
%  further manipulation, e.g. in Illustrator)
%lfp_mtspectrum(..., 'log')
%  Display in decibels, i.e. 10 * log10 of the power spectrum or coherence.
%  This only affects the display, not the value of <Spectra> returned.
%lfp_mtspectrum(..., 'multitrig')
%  Accumulates data over multiple instances of the alignment event per
%  trial.  (The default is to use only the first instance.)  Issues a
%  warning if there are two successive instances that are within
%  lfp_XLimAll(2) - lfp_XLimAll(1) of each other.  Raises an error if
%  lfp_XLimAll and <window> are both empty.
%lfp_mtspectrum(..., 'multiwin', winwidth)
%   Aggregate over umpteen half-overlapped windows of width <winwidth>
%   spanning <window>.  For a given <winwidth>, computation time goes
%   linearly with the duration of <window>, whereas for an equivalently
%   narrow value of 'nw' using the maximum number of tapers, it goes like
%   the square of <window>.  The multiple windows are treated the same as
%   if they were additional sets of trials.
%lfp_mtspectrum(..., 'nodisplay')
%  To control display visibility. By default the figure will be visible.
%  Specifying 'nodisplay' still creates the figure, but with 'Visible'
%  set to 'off'.
%lfp_mtspectrum(..., 'noplot')
%  Does not create a figure window or plot results, returns <hF> = [].
%lfp_mtspectrum(..., 'norefOK')
%  As in lfp_findCommonTime, ignore trials that have no reference event.
%lfp_mtspectrum(..., 'notrunc')
%  As in lfp_findCommonTime, ignore triggers whose inclusion would result
%  in truncation of <interval> to something narrower than lfp_XLimAll.  Has
%  no effect if lfp_XLimAll is empty.
%lfp_mtspectrum(..., 'norefOK'
%   simply ignores any trials that do not contain a reference
%   event; this usually would only be used together with 'mutltitrig'
%   and/or 'evtbounds'.
%lfp_mtspectrum(..., 'norm')
%  Display the power as proportion of total power; title is marked "N"
%  (default "R" for "raw").  If combined with 'avg', normalization is done
%  after averaging, and if 'err' is also specified, the error bars are
%  normalized with respect to total power rather than total error.  If
%  combined with 'log', it is done before taking log, so 0 dB represents
%  the total signal energy.  If combined with 'log', 'avg', and 'err', then
%  the average and error are computed first, then normalized,
%  then converted to log units.  If combined with 'ovr', then each trial is
%  normalized separately.  Does nothing if combined with 'coh'.
%lfp_mtspectrum(..., 'nw', nw)
%  Forces "time-bandwidth product" setting for multitapers to <nw>; default
%  is 1.8.
%lfp_mtspectrum(..., 'p', plevel)
%  Sets the p level for computing the significance threshold in
%  coherence.  Default value is .01.
%lfp_mtspectrum(..., 'pad', N)
%  Unless window width happens already to be 2^k samples long for some k,
%  there is always padding, and N controls how many powers of 2 to skip:
%  N=0 pads to the next greater power of 2, N=1 pads to twice that length,
%  N=2 pads to four times that length, etc.  The default value is N=0.
%lfp_mtspectrum(..., 'phi')
%  Works only with 'coh'.  Displays phase in degrees as a second graph
%  plotted below the first, together with error curves if requested, but
%  NOTE that the error curves for phi are ALWAYS theoretical 95% confidence
%  limits regardless of what error parameters are specified.  Negative
%  phase means that the second channel is lagging the first.
%lfp_mtspectrum(..., 'ovr')
%  Overlays the plots of the specified trials on a single figure.  If 'avg'
%  and 'ovr' are specified together, then the trials are averaged but
%  instead of showing multiple filenums or pairs on multiple subplots, the
%  traces are overlaid on a single plot.
%lfp_mtspectrum(..., 'print')
%  Same as for lfp_disp.
%lfp_mtspectrum(..., 'rmBL', BL)
%  Subtracts the baseline spectrum <BL> (as returned by lfp_BLspectrum)
%  from the power spectrum.  Does not work with 'coh'.  When 
%  scale is 'log', the spectra are converted to dB first, then subtracted.
%  If BL spans a different frequency range from lfp_FreqLim, the overlap is
%  used.  If BL is sampled on a different frequency grid from the
%  spectrogram, then BL is interpolated to match using interp1.
%lfp_mtspectrum(..., 'rmdc')
%  Removes the DC component from wave(s) before computing the coherence or
%  spectrum.
%lfp_mtspectrum(..., 'rmtrend')
%  Detrends the wave(s) before computing the coherence or spectrum.  Since
%  detrending removes the best affine (linear) fit from the data, this
%  operation also entails removing DC, so 'rmdc' is ignored if 'rmtrend' is
%  specified.
%lfp_mtspectrum(..., 'rmEP')
%  Removes the average over trials from each trial in wave(s) before
%  computing the coherence or spectrum.  Only works with 'avg'.
%lfp_mtspectrum(..., 'rotate')
%  Works only with 'arrows'.  Adds 90 degrees to the phase before drawing
%  arrows, so zero phase is straight up instead of straight to the right.
%  Does not alter values plotted using 'phi' option.
%lfp_mtspectrum(..., 'session', sessionname)
%  Supplies a default session name to use when specifying Unique Trial IDs
%  instead of internal trial numbers.
%lfp_mtspectrum(..., 'showtrialnums')
%  Labels figure with trial numbers instead of selection rule
%lfp_mtspectrum(..., 'thresh')
%  Works only with 'coh'.  Plots a horizontal line at the p<.01 statistical
%  significance threshold determined by from Bijan's rule of thumb. 
%lfp_mtspectrum(..., 'xscale', scale), lfp_mtspectrum(..., 'yscale', scale)
%  Works only with 'coh'.  Scale factors for arrows or lines in the x
%  direction and y direction respectively.  Default value is 1.
%lfp_mtspectrum(..., 'unwrap')
%  Works only with 'phi'.  Applies the Matlab 'unwrap' function to the
%  phase before plotting.
%lfp_mtspectrum(..., figbase)
%  figbase is a number to add to the trial number to yield the figure
%  number.
%lfp_mtspectrum(trials)
%  All files are plotted using the default window.  Note that <trials> can
%  be a single integer to plot just one trial (no square brackets needed).
%lfp_mtspectrum
%  Plots all files for all selected trials.  Equivalent to lfp_mtspectrum([]).
%OUTPUTS
%   hF: figure handle
%   Spectra: spectral power or coherence.  If options 'avg' and 'err' are
%       both specified, then <Spectra> is a three-column array where the
%       second and third columns are respectively the lower and upper
%       confidence limits.  Otherwise it is in freqs x plotnum format, with
%       3rd dimension for trials when using 'ovr' or 'avg'. Not affected by
%       'log' or 'rmBL' options (see Spectra2plot).
%   f: frequency points
%   phi: spectral phase to match <Spectra> when 'coh' and 'phi' are
%       specified, in degrees.
%   Spectra2plot: exactly what was last plotted, including the effects of
%       'log' and/or 'rmBL' if specified.

%$Rev: 409 $
%$Date: 2020-04-25 15:13:16 -0400 (Sat, 25 Apr 2020) $
%$Author: dgibson $

global lfp_FileNames lfp_SelectedFiles lfp_TrialStyle lfp_ActiveFilenums ...
    lfp_SamplePeriod lfp_TrialIndex lfp_YPwrLim lfp_XLimAll ...
    lfp_FreqLim lfp_DataDir lfp_AlignmentRef

if nargin <3
    window = [];
elseif ~isempty(window) && ~(isa(window, 'double') ...
        && isequal(size(window), [1 2]) )
    error('lfp_mtspectrum:badwindow', ...
        '<window> must be 1x2 number array.' );
end
if nargin < 2 || isempty(filenums)
    filenums = sort(lfp_ActiveFilenums);
end
if any(~ismember(filenums, lfp_ActiveFilenums))
    error('lfp_mtspectrum:noSuchFilenum', ...
        ['You requested non-existent file numbers: ' ...
            num2str(filenums( ...
            ~ismember(filenums, lfp_ActiveFilenums) )) ]);
end
% Apply lfp_SelectedFiles as a mask to filenums:
filenums = filenums(lfp_SelectedFiles(filenums));

if nargin < 1 || isempty(trials)
    trials = 1:size(lfp_TrialIndex,1);
end

if isempty(window)
    if isempty(lfp_XLimAll)
        window = [0 1];
    else
        window = lfp_XLimAll;
    end
end

phi = [];
session = '';
trialstyle = lfp_TrialStyle;
argnum = 1;
arrowflag = false;
avgflag = false;
axinfo.plotflag = true;
BL = [];
cohflag = false;
commontimeopts = {};
display_str = 'on';
EPflag = false;
errflag = false;
evtbounds = {}; %#ok<NASGU>
figflag = false;
getsamplesopts = {};
K = 1;
lineflag = false;
multitrigflag = false;
normflag = false;
normmark = 'R';
numarrows = 0;
numplots = length(filenums);
numwins = 1;
NW = 1.8;
ovrflag = false;
p = .01;
padfactor = 0;
printflag = false;
phiflag = false;
rmdcflag = false;
rmtrendflag = false;
rmEPflag = false;
rotateflag = false;
threshflag = false;
unwrapflag = false;
vscale = 'lin';
winwidth = diff(window);
xscale = 1;
yscale = 1;
while argnum <= length(varargin)
    if isa(varargin{argnum}, 'double')
        if figflag
            error('lfp_mtspectrum:multifigbase', ...
                'You have specified figbase more than once.');
        else
            figflag = true;
            figbase = varargin{argnum};
        end
    else
        switch varargin{argnum}
            case 'arrows'
                arrowflag = true;
                argnum = argnum + 1;
                numarrows = varargin{argnum};
            case 'lines'
                lineflag = true;
                argnum = argnum + 1;
                numarrows = varargin{argnum};
            case 'avg'
                avgflag = true;
            case 'EP'
                EPflag = true;
            case 'err'
                errflag = true;
                argnum = argnum + 1;
                err = varargin{argnum};
            case 'evtbounds'
                argnum = argnum + 1;
                commontimeopts = [commontimeopts ...
                    {'evtbounds' varargin{argnum}}]; %#ok<*AGROW>
                getsamplesopts = [getsamplesopts ...
                    {'evtbounds' varargin{argnum}}]; %#ok<*AGROW>
            case 'k'
                argnum = argnum + 1;
                K = varargin{argnum};
            case 'log'
                vscale = 'log';
            case 'multitrig'
                multitrigflag = true;
                commontimeopts = [commontimeopts {'multitrig'}];
                getsamplesopts = [getsamplesopts {'multitrig'}];
            case 'multiwin'
                argnum = argnum + 1;
                winwidth = varargin{argnum};
                % Successive windows overlap by one half the window width
                % to accommodate the classic single taper.  Therefore, the
                % number of windows is determined based on half the
                % specified window width.
                numwins = floor(diff(window)/(winwidth/2)) - 1;
            case 'nodisplay'
                display_str = 'off';
            case 'noplot'
                axinfo.plotflag = false;
            case 'notrunc'
                commontimeopts = [commontimeopts {'notrunc' window}];
                getsamplesopts = [getsamplesopts {'notrunc' window}];
            case 'norefOK'
                commontimeopts = [commontimeopts {'norefOK'}];
                getsamplesopts = [getsamplesopts {'norefOK'}];
            case 'norm'
                normflag = true;
                normmark = 'N';
            case 'nw'
                argnum = argnum + 1;
                NW = varargin{argnum};
            case 'coh'
                cohflag = true;
                numplots = length(filenums)/2;
            case 'ovr'
                ovrflag = true;
            case 'p'
                argnum = argnum + 1;
                if ~isa(varargin{argnum}, 'double')
                    error('lfp_mtspectrum:badp', ...
                        'The p level must be a number');
                end
                p = varargin{argnum};
            case 'pad'
                argnum = argnum + 1;
                if ~isa(varargin{argnum}, 'double')
                    error('lfp_mtspectrum:badpad', ...
                        'The padding factor must be a number');
                end
                padfactor = varargin{argnum};
            case 'phi'
                phiflag = true;
            case 'print'
                printflag = true;
            case 'rotate'
               rotateflag = true;
            case 'rmBL'
                argnum = argnum + 1;
                BL = varargin{argnum};
            case 'rmdc'
                rmdcflag = true;
            case 'rmtrend'
                rmtrendflag = true;
            case 'rmEP'
                rmEPflag = true;
            case 'session'
                argnum = argnum + 1;
                session = varargin{argnum};
            case 'showtrialnums'
                trialstyle = 'trialnums';
            case 'thresh'
                threshflag = true;
            case 'unwrap'
                unwrapflag = true;
            case 'xscale'
                argnum = argnum + 1;
                xscale = varargin{argnum};
            case 'yscale'
                argnum = argnum + 1;
                yscale = varargin{argnum};
            otherwise
                error('lfp_mtspectrum:badoption', ...
                    ['The option "' varargin{argnum} ...
                        '" is not recognized.'] );
        end
    end
    argnum = argnum + 1;
end

if ~axinfo.plotflag && nargout < 2
    error('lfp_mtspectrum:noplot', ...
        'You have specified ''noplot'' without return values.');
end

if ischar(trials)
    trials = lfp_parseTrialStr(trials, session);
end
% Apply mask to trials:
trials = lfp_enabledTrials(trials);

if multitrigflag
    if isempty(window) && isempty(lfp_XLimAll)
        error('lfp_mtspectrum:badopt2', ...
            '''multitrig'' requires a pre-set time window');
    end
end

% Ensure that <trials> is an integer row vector
if size(trials,1) > 1
    if size(trials,2) > 1
        error('lfp_mtspectrum:badTrials1', '<trials> must be an integer vector.');
    end
    trials = trials';
end
if ~(isa(trials, 'double')) ...
        || ~all(fix(trials) == trials)
    error('lfp_mtspectrum:badTrials2', '<trials> must be an integer vector.');
end

if ~isequal(class(numarrows), 'double') || numel(numarrows) ~= 1
    error('lfp_mtspectrum:badnumarrows', ...
        'You must specify a number of arrows following the ''arrows'' option');
end

% Check for potential logical errors
if round((window(2) - window(1))/lfp_SamplePeriod) < 2 * NW
    error('lfp_mtspectrum:badwindow2', ...
        'Your window start and stop times %s do not enclose enough samples for nw=%d', ...
        dg_thing2str(window), NW );
end

if any(trials > size(lfp_TrialIndex,1) | trials < 1)
    error('lfp_mtspectrum:badtrialnum', ...
        [ 'You have specified trials that do not exist: ' ...
            num2str(trials( ...
            trials > size(lfp_TrialIndex,1) | trials < 1 )) ]);
end
if cohflag && (mod(length(filenums), 2) ~= 0)
    error('lfp_mtspectrum:badnumfiles', ...
        'You must specify an even number of files for coherence.' );
end
if cohflag && normflag
    normflag = false;
end
if ~cohflag && threshflag
    error('lfp_mtspectrum:badoptcombo3', ...
        '''thresh'' only works with ''coh''');
end
if ~cohflag && arrowflag
    error('lfp_mtspectrum:badoptcombo2', ...
        '''arrows'' only works with ''coh''');
end
if ovrflag && errflag
    error('lfp_mtspectrum:badoptcombo', ...
        'You cannot show error on overlaid traces');
end
if errflag
    if ~ismember(err(1), [ 1 2 ])
        error('lfp_mtspectrum:err', ...
            '<n> value for ''err'' option must be 1 or 2');
    end
    if err(2) > 1 || err(2) < 0
        error('lfp_mtspectrum:err2', ...
            '<p> valuefor ''err'' option must be between 0 and 1 inclusive');
    end
end
if EPflag && avgflag || EPflag && rmEPflag || EPflag && ovrflag
    error('lfp_mtspectrum:EP', ...
        '''EP'' cannot work with ''avg'', ''ovr'', or ''rmEP''.');
end
if avgflag || ovrflag || EPflag
    % Just one figure window
    if figflag
        fignums = figbase;
    else
        fignums = 0;
    end
else    % no trial aggregation
    % A separate figure window for each trial
    if figflag
        fignums = figbase + trials;
    else
        fignums = zeros(size(trials));
    end
end

Fs = 1/lfp_SamplePeriod;

aggrmode = '';
if avgflag
    aggrmode = 'avg';
elseif ovrflag
    aggrmode = 'ovr';
elseif EPflag
    aggrmode = 'EP';
end
axinfo.xlabel = 'Frequency, Hz';
axinfo.xlim = [];
if cohflag && ~isequal(vscale, 'log')
    if ~phiflag
        axinfo.ylim = [0 1];
    else
        for plotnum = 1:numplots
            axinfo.ylim(2*plotnum - 1, :) = [0 1];
            axinfo.ylim(2*plotnum, :) = [0 0];
        end
    end 
else
    axinfo.ylim = lfp_YPwrLim;
end


%
% Done with argument preprocessing.  Actual analysis/display starts here.
%

% Prepare data for chronux funcs in samples x trials form, using 3rd
% dimension for channels
if ~isempty(lfp_XLimAll)
    commontimeopts = [commontimeopts {'recseg'}];
end
[interval, rawtrialinfo] = lfp_findCommonTime(trials, commontimeopts{:});
if isempty(rawtrialinfo)
    warning('lfp_mtspectrum:noevts', ...
        'No alignment events were found.');
    hF = [];
    Spectra = [];
    f = [];
    phi = [];
    ntrigs = [];
    return
end
numtrials = sum(rawtrialinfo(:,3)~=0);
for plotnum = 1:numplots
    if cohflag
        if phiflag
            % Odd-numbered plotnames are coherence, even-numbered are phase
            plotnames{2*plotnum - 1} = ...
                [lfp_FileNames{filenums(2*plotnum - 1)} ...
                ' x ' lfp_FileNames{filenums(2*plotnum)} ];
            plotnames{2*plotnum} = 'phase';
        else
            plotnames{plotnum} = ...
                [lfp_FileNames{filenums(2*plotnum - 1)} ...
                ' x ' lfp_FileNames{filenums(2*plotnum)} ];
        end
    else
        plotnames = lfp_FileNames(filenums);
    end
    % collect arrays of samples to analyze.  Multiple windows created by
    % 'multiwin' option should behave like trials, i.e. end up as
    % additional columns.  We always assign to <data1>, but <data2> is only
    % used rarely, so we initialize it to empty. Due to roundoff errors,
    % the number of samples returned by lfp_getSamples can occasionally
    % differ by 1 sample, so we calculate the start and end of each window
    % <onewin> as an integral number of samples.
    data2 = [];
    winwidthsamples = round(winwidth / lfp_SamplePeriod);
    for winnum = 1:numwins
        onewin = window(1) + lfp_SamplePeriod * ( ...
            floor((winnum - 1) * winwidthsamples / 2) ...
            + [0 winwidthsamples] );
        if cohflag
            data1(:, (winnum-1)*numtrials+1:winnum*numtrials, ...
                plotnum) = lfp_getSamples(trials, ...
                filenums(2*plotnum-1), onewin, getsamplesopts{:});
            data2(:, (winnum-1)*numtrials+1:winnum*numtrials, ...
                plotnum) = lfp_getSamples(trials, ...
                filenums(2*plotnum), onewin, getsamplesopts{:});
        else
            data1(:, (winnum-1)*numtrials+1:winnum*numtrials, ...
                plotnum) = lfp_getSamples(trials, ...
                filenums(plotnum), onewin, getsamplesopts{:});
        end
    end
    % preprocess according to options:
    if rmtrendflag
        data1(:,:,plotnum) = detrend(data1(:,:,plotnum));
        if ~isempty(data2)
            data2(:,:,plotnum) = detrend(data2(:,:,plotnum));
        end
    elseif rmdcflag
        data1(:,:,plotnum) = data1(:,:,plotnum) - ...
            repmat(mean(data1(:,:,plotnum),1),size(data1,1),1);
        if ~isempty(data2)
            data2(:,:,plotnum) = data2(:,:,plotnum) - ...
                repmat(mean(data2(:,:,plotnum),1),size(data2,1),1);
        end
    end
    if rmEPflag
        data1(:,:,plotnum) = data1(:,:,plotnum) - ...
            repmat(mean(data1(:,:,plotnum),2),1,size(data1,2));
        if ~isempty(data2)
            data2(:,:,plotnum) = data2(:,:,plotnum) - ...
                repmat(mean(data2(:,:,plotnum),2),1,size(data2,2));
        end
    end
    if isequal(aggrmode, 'EP')
        if plotnum == 1
            meandata1 = NaN(size(data1, 1), size(data1, 3));
            if ~isempty(data2)
                meandata2 = NaN(size(data2, 1), size(data2, 3));
            end
        end
        meandata1(:,plotnum) = mean(data1(:,:,plotnum),2);
        if ~isempty(data2)
            meandata2(:,plotnum) = mean(data2(:,:,plotnum),2);
        end
    end
end
ntrigs = size(data1, 2);
if isequal(aggrmode, 'EP')
    data1 = meandata1;
    if ~isempty(data2)
        data2 = meandata2;
    end
end

if cohflag
    % The following formula is based on the approximation that the
    % distribution of coherence is Gaussian (Jarvis & Mitra 2001 section
    % 5.1) centered at zero (it is actually more like 0.05, from Jarvis &
    % Mitra 2001 eqn 5.4, but we normally look at coherence >> .05 so we
    % call it 0) with variance 1/sqrt(nu) where nu is the number of
    % independent observations (Mitra and Pesaran 1999).
    thresh = tinv(1-p,K*ntrigs)/sqrt(K*ntrigs);
end

% Save results into Spectra as freqs x plotnum (3rd dimension is trials for
% 'ovr' or 'avg' option).
for plotnum = 1:numplots
    if isempty(lfp_FreqLim)
        freqlim = [0 1/(2*lfp_SamplePeriod)];
    else
        freqlim = lfp_FreqLim;
    end
    if ~isempty(BL)
        if ~isequal(BL.f([1 end]), freqlim)
            freqlim(1) = max(freqlim(1), BL.f(1));
            freqlim(2) = min(freqlim(2), BL.f(end));
        end
    end
    if errflag
        if cohflag
            if err(1)==2
                [Spectra(:,plotnum,:),phi(:,plotnum,:),f,confC, ...
                    phierr(:,plotnum,:),...
                    errcurve(1:2,:,plotnum,:)] = chronux_coherencyc( ...
                    data1(:,:,plotnum), data2(:,:,plotnum), ...
                    [NW K], padfactor, Fs, freqlim, err, avgflag); %#ok<*ASGLU>
            else
                % show stat. signif. threshold:
                [Spectra(:,plotnum,:),phi(:,plotnum,:),f,confC, ...
                    phierr(:,plotnum,:)] = chronux_coherencyc( ...
                    data1(:,:,plotnum), data2(:,:,plotnum), ...
                    [NW K], padfactor, Fs, freqlim, err, avgflag);
                errcurve(1,:,plotnum,:) = ...
                    ones(size(Spectra(:,plotnum,:))) * confC(1);
                errcurve(2,:,plotnum,:) =  errcurve(1,:,plotnum,:);
            end
        else
            [Spectra(:,plotnum,:),f,errcurve(1:2,:,plotnum,:)] = ...
                chronux_mtspectrumc( data1(:,:,plotnum), ...
                [NW K], padfactor, Fs, freqlim, err, avgflag);
        end
    else    % no error measures, might have 'ovr'
        if cohflag
            [Spectra(:,plotnum,:), phi(:,plotnum,:), f] = chronux_coherencyc( ...
                data1(:,:,plotnum), data2(:,:,plotnum), ...
                [NW K], padfactor, Fs, freqlim, 0, avgflag);
        else
            [Spectra(:,plotnum,:), f] = chronux_mtspectrumc( ...
                data1(:,:,plotnum), ...
                [NW K], padfactor, Fs, freqlim, 0, avgflag);
        end
    end
    if ~isempty(BL)
        BLidx = (BL.f >= freqlim(1) & BL.f <= freqlim(2));
        if length(BL.f(BLidx)) ~= length(f)
            intrpBL = reshape( interp1( BL.f(BLidx), ...
                BL.sum(BLidx), f ), [], 1) / BL.N ;
        else
            intrpBL = BL.sum(BLidx)/BL.N;
        end
    end
end

% Produce figs
if ovrflag && ~avgflag
    option = 'ovr';
elseif ovrflag && avgflag
    option = 'ovrchan';
elseif errflag
    option = 'asymmetric';
else
    option = '';
end
for figwinnum = 1:length(fignums)
    if avgflag || ovrflag || EPflag
        trialslabel = lfp_getTrialsLabel(trials, trialstyle);
        % If cohflag, then Spectra can only be 3D here if ovrflag (note
        % that this does not execute if ~(avgflag || ovrflag || EPflag)).
        if avgflag || EPflag
            Spectra2plot = mean(Spectra, 3);
        else    % ovrflag & ~avgflag
        % Spectra2plot 3rd dimension: trials 
            Spectra2plot = Spectra;
        end
        if cohflag
            % This is NOT where phi gets aggregated!
            Phases2plot = phi;
        end
    else    % no aggregation, so 1 trial per fig
        trialslabel = lfp_getTrialsLabel(trials(figwinnum), trialstyle);
        Spectra2plot = Spectra(:,:,figwinnum);
        if cohflag
            Phases2plot = phi(:,:,figwinnum);
        end
    end
    trialslabel = ['trials ' trialslabel];
    if multitrigflag
        ntrigstr = sprintf(' ntrigs=%d', ntrigs);
    else
        ntrigstr = '';
    end
    if rmtrendflag
        rmstr = 'rmtrend';
    elseif rmdcflag
        rmstr = 'rmdc';
    else
        rmstr = '';
    end
    if rmEPflag
        if isempty(rmstr)
            rmstr = 'rmEP';
        else
            rmstr = [rmstr ' rmEP'];
        end
    end
    figtitle = sprintf('%s align=%s win=%s %s %s %s\n%s %s%s\nnw=%g k=%g pad=%g', ...
        lfp_DataDir, mat2str(lfp_AlignmentRef), mat2str(window), ...
        rmstr, vscale, normmark, aggrmode, ...
        trialslabel, ntrigstr, ...
        NW, K, padfactor);
    if numwins > 1
        figtitle = sprintf('%s; multiwin:%g', figtitle, winwidth);
    end
    if cohflag
        figtitle = sprintf('%s; p=%g level: %g', figtitle, p, thresh);
    end
    if unwrapflag
        Phases2plot(:,:,1) = unwrap(Phases2plot(:,:,1));
    end
    if errflag
        % Spectra2plot 3rd dimension:
        % 1=mean 2=lowerbound 3=upperbound 
        Spectra2plot(:,:,2) = errcurve(1,:,:,figwinnum);
        Spectra2plot(:,:,3) = errcurve(2,:,:,figwinnum);
        if cohflag
            % see coherr.m and Jarvis & Mitra section 5.2 regarding the
            % following calculation of 95% confidence limits:
            Phases2plot(:,:,2) = Phases2plot(:,:,1) - 2 * phierr(:,:,figwinnum);
            Phases2plot(:,:,3) = Phases2plot(:,:,1) + 2 * phierr(:,:,figwinnum);
            figtitle = [ figtitle sprintf(' err=%s', dg_thing2str(err)) ];
        end
    end
    if normflag
        if ovrflag
            for trialidx = 1:size(Spectra2plot,3)
                for plotnum = 1:size(Spectra2plot,2)
                    Spectra2plot(:,plotnum,trialidx) = ...
                        Spectra2plot(:,plotnum,trialidx) ...
                        / sum(Spectra2plot(:,plotnum,trialidx));
                end
            end
        else
            for plotnum = 1:size(Spectra2plot,2)
                Spectra2plot(:,plotnum,:) = ...
                    Spectra2plot(:,plotnum,:) ...
                    / sum(Spectra2plot(:,plotnum,1));
            end
        end
    end
    if isequal(vscale, 'log')
        Spectra2plot = 10 * log10(Spectra2plot);
        if ~isempty(BL)
            intrpBL = 10*log10(intrpBL);
        end
    end
    if ~isempty(BL)
        Spectra2plot = Spectra2plot - repmat(intrpBL, ...
            [1 size(Spectra2plot,2) size(Spectra2plot,3)] );
    end
    % Spectra2plot is freqs x plotnum x (see comments above).
    % At this point, we finally interleave coherence and phase plots if
    % necessary:
    if cohflag && phiflag
        for k = 1:size(Spectra2plot,2)
            curves2plot(:, 2*k-1, :) = Spectra2plot(:, k, :);
            curves2plot(:, 2*k, :) = 360*Phases2plot(:, k, :) / (2*pi);
        end
    else
        curves2plot = Spectra2plot;
    end
    if axinfo.plotflag
        if fignums(figwinnum) == 0
            fignums(figwinnum) = figure('Visible', display_str);
        else
            set(fignums(figwinnum), 'Visible', display_str);
        end
        hF = lfp_multichannel_plot(fignums(figwinnum), ...
            figtitle, plotnames, ...
            curves2plot, ...
            f', 0, axinfo, option);
        hA = findobj(hF, 'Type', 'axes');
        for k = 1:length(hA)
            % <hA> is in reverse order, so flip it to make <graphnum>
            % start at top:
            graphnum = length(hA) - k + 1;
            if ~phiflag || mod(graphnum, 2)
                % <plotn> is meant to correspond to <plonum> in the earlier
                % loop:
                if phiflag
                    plotn = (graphnum - 1)/2 + 1;
                else
                    plotn = graphnum;
                end
                % Put extra markings on coherence plots, not phase plots:
                if threshflag
                    xlimits = get(hA(k), 'XLim');
                    hold(hA(k), 'on');
                    plot(hA(k), xlimits, [thresh thresh],'r');
                end
                if arrowflag || lineflag
                    xlimits = get(hA(k), 'XLim');
                    ylimits = get(hA(k), 'YLim');
                    xsize = (xlimits(2) - xlimits(1)) / numarrows;
                    ysize = (ylimits(2) - ylimits(1)) / numarrows;
                    freqs = ((1:numarrows)' - 0.5) * xsize;
                    cohvals = interp1(f', Spectra2plot(:,plotn), freqs);
                    phivals = interp1(f', phi(:,plotn,1), freqs);
                    signifmask = cohvals > thresh;
                    if rotateflag
                        phivals = phivals + pi/2;
                    end
                    if any(signifmask)
                        clear starts stops;
                        starts(:,1) = freqs(signifmask);
                        starts(:,2) = cohvals(signifmask);
                        stops(:,1) = starts(:,1) + xscale*xsize*cos(phivals(signifmask));
                        stops(:,2) = starts(:,2) + yscale*ysize*sin(phivals(signifmask));
                        hold(hA(k), 'on');
                        set(0, 'CurrentFigure', hF);
                        set(hF, 'CurrentAxes', hA(k));
                        if arrowflag
                            % Note: Benedikt Halldorsson reported problems with certain
                            % linestyles when printing these arrows, so beware
                            % (Erik Johnson's function downloaded from
                            % http://www.mathworks.com/matlabcentral/fileexchange )
                            arrow(starts, stops, 'Length', round(16*(10/numarrows)));
                        else    % must be lineflag
                            for linenum = 1:size(starts,1)
                                plot([starts(linenum,1) stops(linenum,1)], ...
                                    [starts(linenum,2) stops(linenum,2)], 'g' );
                            end
                        end
                    end
                end
            end
        end
        if printflag
            print;
        end
    else
        hF = [];
    end
end

if errflag && avgflag
    Spectra = [Spectra(:,1,1) errcurve(1:2,:,1,1)'];
    if cohflag && phiflag
        phi = [curves2plot(:,2,1) curves2plot(:,2,2) curves2plot(:,2,3)];
    end
end

