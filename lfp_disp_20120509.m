function [hF, data, ntrigs] = lfp_disp(trials, filenums, window, varargin)
% euh... maybe not:
% lfp_disp:
%     changed behavior of 'norm' option when used with 'xcov'
%     added code to skip unnecessary computation of non-displayed lags when
%         using 'avgxcov' together with non-empty lfp_Xxcov
% 09-May-2012 17:48:59 abandoned development of this in favor of adding
% lfp_getSamples, lfp_xcov, etc.
        
        %modified from v.85 by JF & TF
%LFP_DISP shows waveforms or their crosscovariance functions.
%lfp_disp(trials, filenums, window)
%lfp_disp(..., 'avg')
%lfp_disp(...,'avg','rect')
%lfp_disp(..., 'avg', 'err')
%lfp_disp(..., 'avg', 'err2')
%lfp_disp(..., 'avg', 'err3')
%lfp_disp(..., 'avg', 'err4', numboots)
%lfp_disp(..., 'avg', 'file')
%lfp_disp(..., 'avg', 'fileflip')
%lfp_disp(..., 'avg', 'dB', dBref)
%lfp_disp(..., 'bigevts', bigevtcodes)
%lfp_disp(..., 'evtavg', evts2avg)
%lfp_disp(..., 'evtavg2', evts2avg)
%lfp_disp(..., 'evtbounds', {startIDs stopIDs})
%lfp_disp(..., 'filenames', {<name> ...})
%lfp_disp(..., 'avgxcov')
%lfp_disp(..., 'avgxcov', 'shuffle')
%lfp_disp(..., 'eye')
%lfp_disp(..., 'eye', 'cal')
%lfp_disp(..., figbase)
%lfp_disp(..., 'marknanbad')
%lfp_disp(..., 'markrefOK')
%lfp_disp(..., 'mov')
%lfp_disp(..., 'multi')
%lfp_disp(..., 'multising')
%lfp_disp(..., 'multitrig')
%lfp_disp(..., 'multiwindow')
%lfp_disp(..., 'nodisplay')
%lfp_disp(..., 'noplot')
%lfp_disp(..., 'norefOK')
%lfp_disp(..., 'ovr')
%lfp_disp(..., 'pass')
%lfp_disp(..., 'print')
%lfp_disp(..., 'save')
%lfp_disp(..., 'rmEP')
%lfp_disp(..., 'xcov')
%lfp_disp(..., 'xcov', 'norm')
%lfp_disp(..., 'session', sessionname)
%lfp_disp(..., 'showtrialnums')
%lfp_disp(trials, filenums)
%lfp_disp(trials)
%lfp_disp

%lfp_disp(trials, filenums, window)
%  Plots CSC data from the file numbers in filenums together with all
%  events that fall within the trial boundaries (i.e. between
%  lfp_NominalTrialStart and lfp_NominalTrialEnd).  Each enabled trial
%  (i.e. lfp_SelectedTrials(trial) is 'true' and <trial> is not a member of
%  lfp_BadTrials) in the integer vector <trials> is plotted in its own
%  figure window, numbered with the trial number.  <trials> can be a row
%  vector or a column vector, or it can be empty, in which case all
%  selected trials are used (equivalent to specifying <trials> as
%  find(lfp_SelectedTrials)). If <trials> is a string, then it is
%  interpreted as containing Unique Trial IDs, i.e. the combination of a
%  session name and a trial number as it was originally numbered in that
%  session.  The session name component is optional (see lfp_parseTrialStr
%  and lfp_getTrialNum for further details). Time is shown relative to the
%  first lfp_AlignmentRef event between the trial's lfp_NominalTrialStart
%  and lfp_NominalTrialEnd events.  If there is no lfp_AlignmentRef event,
%  then absolute time is shown.  If <filenums> is the empty array, then all
%  files are plotted. If <window> is specified and not empty, it overrides
%  lfp_XLimAll.  If it is a 1x2 array of double, then the value is used
%  directly.  If it is a 1x2 cell array, then the first cell should contain
%  a list of event IDs, and the event in the list that occurs ealiest in
%  the trial is used as the start time.  The second cell works similarly to
%  specify end time.  Note: You cannot specify 'avg' or 'ovr' with an
%  event-bounded window.
%lfp_disp(..., figbase)
%  The figure window number is assigned to be <figbase> + <trial>.
%lfp_disp(..., 'avg')
%lfp_disp(..., 'avg', 'err')
%lfp_disp(..., 'avg', 'err2')
%lfp_disp(..., 'avg', 'err2', 'dB', dBref)
%lfp_disp(..., 'avg', 'err3', signif)
%lfp_disp(..., 'avg', 'err4', numboots)
%  Instead of plotting each trial in its own figure window, the data for
%  the trials are averaged first and then plotted.  In this case, the only
%  event marker shown is for the lfp_AlignmentRef event.  We make the
%  approximation here that lfp_SamplePeriod is exact and constant,
%  regardless of how many frames of CSC data are included in the time range
%  to be averaged.  Adding 'err' shows additional curves at avg + sd and
%  avg - sd, where sd is the standard deviation at each time point.  'err2'
%  is the same but shows twice the standard error of the mean (sd/sqrt(N))
%  instead of standard deviation.  'err3' shows curves at the positive and
%  negative confidence limits for zero-average noise at the p < <signif>
%  level determined using Student's t distribution. 'err4' shows curves at
%  the positive and negative confidence limits computed by the bootstrap
%  method with <numboots> iterations.  This can become quite
%  computation-intensive as <numboots> increases, and will usually not
%  yield substantially different results from 'err2', owing to the fact
%  that the average of even a few trials will be distributed approximately
%  normally even when the individual trials are not.  'dB' only works in
%  conjunction with 'err2', and does not affect the calculation in any way,
%  but the final plots are converted to dB using 10*log10*(result) -
%  <dbref> before plotting.
%lfp_disp(..., 'bigevts', bigevtcodes)
%  As in lfp_spikeAnalysis.
%lfp_disp(..., 'evtavg', evts2avg)
%  Only works with 'avg' or 'ovr'.  Collects relative event times of all
%  events during the trial whose events whose IDs are in <evts2avg> and
%  plots an event marker at the median event time.  The clickable info for
%  the event marker includes median, mean, SD, first quartile, last
%  quartile, 5th percentile, 95th percentile.  Due to roundoff errors when
%  converting between event times and sample numbers, temporal resolution
%  is limited to half a sample period (lfp_SamplePeriod/2) on each trial,
%  although resolution improves when averaging over many trials; it is
%  unpredictable as to whether the lfp_NominalTrialStart and
%  lfp_NominalTrialEnd events will be included as part of the trial or not
%  (this depends on whether they round to the next higher or next lower
%  sample number); and there is a bug that may cause erratic behavior if
%  lfp_NominalTrialStart is used as the value of lfp_AlignmentRef.
%lfp_disp(..., 'evtavg2', evts2avg)
%  Same as 'evtavg', except only events that fall within the trial
%  boundaries AND the display window are used.
%lfp_disp(..., 'evtbounds', {startIDs stopIDs})
%  Sets the range of time to search for alignment events to something
%  narrower than the entire trial, specified in terms of other events.
%  <startIDs> is a list of alternative event IDs to use as the start of the
%  time range, and <endIDs> is a list of alternative event IDs to use as
%  the end of the time range.  The start event is the earliest event in the
%  trial that is a member of <startIDs>; the end event is the earliest
%  event in the trial AFTER the start event that is a member of <endIDs>.
%  This accomplishes the same effect as giving {startIDs stopIDs} as the
%  value of the <window> parameter.  Note: as of 7/25/05, this does not
%  work with Joey eye plots (probably Theresa too).
%lfp_disp(..., 'eye')
%lfp_disp(..., 'eye', 'cal')
%  Show 2-D plot of eye position with time shown as color.  <filenums> must
%  have either two or three elements.  In either case, the first channel is
%  used as X position and the second as Y position. If there is no third
%  channel given, then eye traces are plotted using uniform size dots.  If
%  there is a third channel, then lfp_plot_eye_t makes dot diameter
%  proportional to the value of the signal on the third channel.  If 'cal'
%  is specified, the eye traces are assumed to be already calibrated;
%  default assumes they are raw eye tracker files.  This will also work
%  with video tracker data from rodent experiments (works with
%  lfp_SetupType = 'monkey', 'theresa', 'rodent').  Behavior is undefined
%  when combined with 'ovr' or 'avg' (for lfp_SetupType = 'rodent',
%  specifying multiple trials to plot automatically overlays).
%lfp_disp(..., 'eye', 'trackeronly')
%  Suppresses setup-specific plotting and just plots the tracker data as
%  black points.
%lfp_disp(..., 'file', N)
%lfp_disp(..., 'fileflip', N)
%  Adding 'file' to 'avg' causes lfp_disp to save the tabulated per-trial
%  data to a spreadsheet at the same point in processing where the 'err'
%  values are calculated.  ('file' can be combined with 'err', too.)  Note
%  that a separate spreadsheet must be saved for each filenum on display.
%  'fileflip' is the same as 'file', except that the rows and columns are
%  flipped, i.e. each column shows one trial.  This is necessary whenever
%  the length of data from one trial exceeds 256 data points.  N is an
%  integer that specifies how densely to sample the data points, i.e. one
%  out of every N data points will be written to the file, starting with
%  the first data point (no data points are skipped at the beginning).
%  NOTES: this does not work with 'xcov', but instead puts the original raw
%  waveform data in the spreadsheet; also, although it works with
%  'avgxcov', only the first half of the xcov function shows up in
%  spreadsheet.
%lfp_disp(..., 'filenames', {<name> ...})
%  For options that save data to files (e.g. 'file'), specifies the
%  filenames to save to.  Must be followed by a cell string array of
%  filenames.  If this option is not given, then a GUI is presented for
%  entering filenames.
%lfp_disp(...,'avg','rect')
%  Rectifies the waveforms before averaging them.
%lfp_disp(..., 'avgxcov')
%  Computes xcovs first, then averages them.  Can be combined with 'err' to
%  get +/- sd curves, and with 'norm' to get xcov(... 'coeff').
%lfp_disp(..., 'avgxcov', 'shuffle')
%  Computes the cross-covariance between two channels, as for 'avgxcov',
%  but the trials of the second channel are randomly shuffled s.t. no trial
%  is in its original position (cf. the 'shuffle' option in lfp_spec)
%lfp_disp(..., 'marknanbad')
%  When called with 'avg', adds to lfp_BadTrials any trials that contain
%  NaN samples in the time window being averaged.  May crash if NOT called
%  with 'avg'.  If ANY channel in <filenums> contains NaN, the trial is
%  marked bad.
%lfp_disp(..., 'markrefOK')
%  Same as 'norefOK', except that it marks the trials that did have
%  reference events by setting lfp_SelectedTrials 'true' for those trials
%  and 'false' for all other trials.  WARNING:  USE OF THIS OPTION ALTERS
%  THE VALUE OF THE GLOBAL <lfp_SelectedTrials>.
%lfp_disp(..., 'nodisplay')
%  To control display visibility. By default the figure will be visible.
%  Specifying 'nodisplay' still creates the figure, but with 'Visible'
%  set to 'off'.
%lfp_disp(..., 'multi', [m_plots_t n_plots_t])
%  For use only with calls to lfp_plot_eye_t, plots <m_plots_t>*<n_plots_t>
%  in a single figure.
%lfp_disp(..., 'multising', [m_plots n_plots])
%  Plots <m_plots>*<n_plots> trials for a single file in a single figure
%lfp_disp(..., 'multitrig')
%  When combined with 'avg', enables event-triggered averaging over
%  multiple occurences of the event per trial.  Bombs out if lfp_XLimAll
%  and <window> are both empty. Appears to work for 'ovr' as well, but I'm
%  not sure this has been tested.
%lfp_disp(..., 'multiwindow')
%  Displays 4 different time windows per trial in separate subplots.  Only
%  works with Joey 2D eye plots (and possibly Theresa).
%lfp_disp(..., 'noplot')
%  If 'avg' or 'ovr' is specified, skips plotting and returns [] for <hF>.
%  Has no effect if neither 'avg' nor 'ovr' is specified.
%lfp_disp(..., 'norefOK')
%  Simply skips any trials trials that do not have a reference event.
%lfp_disp(..., 'ovr')
%  Overlays traces on a single figure.  If 'avg' and 'ovr' are specified
%  together, then the trials are averaged but instead of showing multiple
%  filenums or pairs on multiple subplots, the traces are overlaid on a
%  single plot ('ovrchan' plot); in this case if there is only one trial to
%  plot, then event markers are drawn as when 'avg' is not given.  Clears
%  lfp_ClickedTrials.
%lfp_disp(..., 'pass', {arguments})
%  use this option to pass arguments on to 'eye' plotting functions for
%  lfp_SetupType in {'rodent' 'theresa'} (i.e. lfp_plot_eye_rodent and
%  lfp_plot_eye_t).
%lfp_disp(..., 'print')
%  As soon as the figure has finished rendering, it is printed on the
%  default printer.  If you change your default printer, then you must
%  restart Matlab in order for the change to take effect.
%lfp_disp(..., 'rmEP')
%  Only works when used together with 'avgxcov' or 'ovr'.  Once the array
%  of samples from all the individual trials has been constructed, it is
%  averaged across trials and the trial average is subtracted from each
%  individual trial.
%lfp_disp(..., 'save')
%  Automatically saves as a *.fig file to an automatically generated
%  filename.
%lfp_disp(..., 'session', sessionname)
%  Supplies a default session name to use when specifying Unique Trial IDs
%  instead of internal trial numbers.
%lfp_disp(..., 'showtrialnums')
%  Labels figure with trial numbers instead of selection rule
%lfp_disp(..., 'xcov')
%lfp_disp(..., 'xcov', 'norm')
%  Computes the cross-covariance between two channels. <filenums> must
%  have an even number of entries after being masked by lfp_SelectedFiles;
%  a separate plot will be shown for each pair. If 'norm' is specified,
%  then dg_xcov(x, y, 'coeff') is called, if not then xcov(x, y, 'unbiased').
%  When 'xcov' is combined with 'avg', the waveforms are averaged first,
%  and then submitted to xcov.  <window> and lfp_XLimAll are applied to the
%  waveforms first before submitting to xcov, so what you see without
%  'xcov' is what gets used with 'xcov'.  A large peak at a negative lag in
%  the xcov implies that the first waveform is leading the second.
%lfp_disp(trials)
%  All files are plotted.  Note that <trials> can be a single integer to
%  plot just one trial (no square brackets needed).
%lfp_disp
%  Plots all files for all selected trials.  Equivalent to
%  lfp_disp([]).
%
% OUTPUTS
%[hF, data, ntrigs] = lfp_disp(...)
%   hF: for strip charts, contains hF returned by lfp_multichannel_plot;
%       otherwise [].
%   data: for 'ovr' strip charts, returns the data points in samples X
%       channels X trials form; for other strip charts, contains data
%       returned by lfp_multichannel_plot (see below); for eye plots,
%       contains whatever is returned by the lfp_plot_eye_* function;
%       otherwise [].  What lfp_multichannel_plot returns for certain
%       popular options is:
%           'avg' - a cell array containing a 2-column array of (X,Y)
%               coordinates for each graph plotted.
%           'avg', 'err' - like 'avg' but a 3rd column is added for the
%               error, so conf limits would be data{k}(:,2)+data{k}(:,3)
%               and  data{k}(:,2)-data{k}(:,3).
%           'avg', 'err2' - same as 'avg', 'err'.
%   ntrigs: if multitrig, the actual number of trigger events found;
%       otherwise 0.

%$Rev: 269 $
%$Date: 2012-05-07 15:07:49 -0400 (Mon, 07 May 2012) $
%$Author: dgibson $

lfp_declareGlobals;

colorflipped = false;
hF = [];
data = [];
ntrigs = [];

evtbounds = {};
if nargin <3 || isempty(window)
    window = [];
else
    switch class(window)
        case 'double'
            if ~isequal(size(window), [1 2])
                error('lfp_disp:badwindow1', ...
                    '<window> must be 1x2 number or cell array.' );
            end
        case 'cell'
            if ~isequal(size(window), [1 2])
                error('lfp_disp:badwindow2', ...
                    '<window> must be 1x2 number or cell array.' );
            end
            evtbounds = window;
            window = [];
        otherwise
            if ~isequal(size(window), [1 2])
                error('lfp_disp:badwindow1', ...
                    '<window> must be a number or cell array.' );
            end
    end
end

if nargin < 2 || isempty(filenums)
    filenums = sort(lfp_ActiveFilenums);
end
% Apply lfp_SelectedFiles as a mask to filenums:
filenums = filenums(find(lfp_SelectedFiles(filenums)));

if nargin < 1 || isempty(trials)
    trials = 1:size(lfp_TrialIndex,1);
end
session = '';
trialstyle = lfp_TrialStyle;
argnum = 1;
avgflag = false;
avgxcovflag = false;
shuffleflag = false;
bigevts = {};
calflag = false;
errflag = false;
err2flag = false;
err3flag = false;
err4flag = false;
erropt = '';
evtavg2flag = false;
evts2avg = [];
eyeflag = false;
figflag = false;
fileflag = false;
filenames = {};
flipflag = false;
histflag = false;
joyflag = false;
markrefOKflag = false;
marknanbadflag = false;
multitrigflag = false;
multiwindowflag = false;
normflag = false;
norefOKflag = false;
ovrflag = false;
plotflag = true;
printflag = false;
rectflag = false;
trigfuncflag = false;
xcovflag = false;
saveflag = false;
dispflag = false;
multiflag = false;
passflag = false;
argstoeye = {};
passoverlayflag = false;
rmEPflag = false;
movflag = false;
multisingflag = false; %TMD 2/18/08
notruncflag = false;
display_str = 'on';
trackeronlyflag = false;
trialinfo = [];
numboots = [];
dBflag = false;
display_str = 'on';
refisOK = false(size(lfp_SelectedTrials));

while argnum <= length(varargin)
    if strcmp(class(varargin{argnum}), 'double')
        if figflag
            error('lfp_disp:multifigbase', ...
                'You have specified figbase more than once.');
        else
            figflag = true;
            figbase = varargin{argnum};
        end
    else
        switch varargin{argnum}
            case 'avg'
                avgflag = true;
                trialinfo = [];
            case 'avgxcov'
                avgflag = true;
                trialinfo = [];
                xcovflag = true;
                numplots = length(filenums)/2;
                avgxcovflag = true;
            case 'shuffle'
                shuffleflag = true;
            case 'bigevts'
                argnum = argnum + 1;
                bigevts = varargin{argnum};
                try
                    bigevtIDs = cell2mat(bigevts(:,1));
                catch
                    error('lfp_spikeAnalysis:badbigevts', ...
                        'Value for <bigevts> is badly formatted.' );
                end
            case 'evtavg'
                argnum = argnum + 1;
                evts2avg = varargin{argnum};
            case 'evtavg2'
                argnum = argnum + 1;
                evts2avg = varargin{argnum};
                evtavg2flag = true;
            case 'eye'
                eyeflag = true;
                multitracestartsamples = [];
                multitraceendsamples = [];
            case 'cal'
                calflag = true;
            case 'err'
                errflag = true;
                erropt = varargin{argnum};
            case 'err2'
                err2flag = true;
                erropt = varargin{argnum};
            case 'err3'
                err3flag = true;
                erropt = varargin{argnum};
                error('''err3'' not implemented yet');
            case 'err4'
                err4flag = true;
                erropt = varargin{argnum};
                argnum = argnum + 1;
                numboots = varargin{argnum};
            case 'evtbounds'
                argnum = argnum + 1;
                evtbounds = varargin{argnum};
            case 'file'
                avgflag = true;
                trialinfo = [];
                fileflag = true;
                argnum = argnum + 1;
                filesampling = varargin{argnum};
            case 'fileflip'
                avgflag = true;
                trialinfo = [];
                fileflag = true;
                flipflag = true;
                argnum = argnum + 1;
                filesampling = varargin{argnum};
            case 'filenames'
                argnum = argnum + 1;
                filenames = varargin{argnum};
            case 'joy'
                joyflag = true;
            case 'marknanbad'
                marknanbadflag = true;
            case 'markrefOK'
                markrefOKflag = true;
                norefOKflag = true;
            case 'mov'
                movflag = true;
            case 'multitrig'
                multitrigflag = true;
            case 'multiwindow'
                multiwindowflag = true;
            case 'norefOK'
                norefOKflag = true;
            case 'notrunc'
                notruncflag = true;
            case 'ovr'
                ovrflag = true;
                trialinfo = [];
            case 'noplot'
                plotflag = false;
            case 'rect'
                rectflag = true;
            case 'rect'
                rectflag = true;
            case 'rmEP'
                rmEPflag = true;
            case 'showtrialnums'
                trialstyle = 'trialnums';
            case 'trackeronly'
                trackeronlyflag = true;
            case 'xcov'
                xcovflag = true;
                numplots = length(filenums)/2;
            case 'norm'
                normflag = true;
            case 'print'
                printflag = true;
            case 'session'
                argnum = argnum + 1;
                session = varargin{argnum};
            case 'save'
                saveflag = true;
            case 'disp'
                dispflag = true;
            case 'multi'
                multiflag = true;
                argnum = argnum + 1;
                m_plots_t = varargin{argnum}(1);
                n_plots_t = varargin{argnum}(2);
            case 'pass'
                % use this option to pass arguments on to eye plotters
                passflag = true;
                argnum = argnum + 1;
                argstoeye = varargin{argnum};
                % later, need to know if 'overlay' was one of the options
                % passed so search for it and set flag here:
                arg = 1;
                while arg <= size(argstoeye,2) && ~passoverlayflag
                    if isequal(argstoeye{arg},'overlay')
                        overlay_idx = arg;
                        passoverlayflag = true;
                    end
                    arg = arg + 1;
                end
            case 'multising' %TMD 4/4/08
                if length(filenums) > 1
                    error('lfp_disp:multising', ...
                        'This option should only have one file number.');
                end
                multisingflag = true;
                m_plot_num = 1;
                argnum = argnum + 1;
                m_plots = varargin{argnum}(1);
                n_plots = varargin{argnum}(2);
            case 'dB'
                dBflag = true;
                argnum = argnum + 1;
                dBref = varargin{argnum};
            case 'nodisplay'
                display_str = 'off';
            otherwise
                error('lfp_disp:badoption', ...
                    ['The option "' varargin{argnum} '" is not recognized.'] );
        end
    end
    argnum = argnum + 1;
end

if strcmp(class(trials), 'char')
    trials = lfp_parseTrialStr(trials, session);
end
% <trials> is now numeric, i.e. trialnums.
if any(trials > length(lfp_SelectedTrials))
    warning('lfp_disp:trials', ...
        'Ignoring trials %s, which are beyond the last trial.', ...
        dg_canonicalSeries(trials(trials > length(lfp_SelectedTrials))) );
    trials(trials > length(lfp_SelectedTrials)) = [];
end
% Apply mask to trials:
trials = lfp_enabledTrials(trials);

if any(~ismember(filenums, lfp_ActiveFilenums))
    error('lfp_disp:noSuchFilenum', ...
        ['You requested non-existent file numbers: ' ...
        num2str(filenums(find( ...
        ~ismember(filenums, lfp_ActiveFilenums) ))) ]);
end
if size(trials,1) > 1
    if size(trials,2) > 1
        error('lfp_disp:badTrials1', '<trials> must be an integer vector.');
    end
    trials = trials';
end
if any(trials > size(lfp_TrialIndex,1) | trials < 1)
    error('lfp_disp:badtrialnum', ...
        [ 'You have specified trials that do not exist: ' ...
        num2str(trials(find( ...
        trials > size(lfp_TrialIndex,1) | trials < 1 ))) ]);
end
if ~(strcmp(class(trials), 'double')) ...
        || ~all(fix(trials) == trials)
    error('lfp_disp:badTrials2', '<trials> must be an integer vector.');
end

if fileflag && (~strcmp(class(filesampling), 'double'))
    error(lfp_disp:badfilesampling', ...
        'The ''file'' and ''fileflip'' options must be followed by a number.' );
end

if xcovflag && (mod(length(filenums), 2) ~= 0)
    error('lfp_disp:badnumfiles', ...
        'You must specify an even number of files for crosscovariance.' );
end

if eyeflag && ~ismember(length(filenums), [2 3])
    error('lfp_disp:badnumfiles', ...
        'You must specify two (or optionally three) files for eye plots.' );
end

if errflag + err2flag + err3flag > 1
    error('lfp_disp:multierrs', ...
        'You must choose just one of ''err'', ''err2'', ''err3''.' );
end

if ~isempty(filenames) && (...
        ~isequal(length(filenames), length(filenums)) ...
        || ~isequal(class(filenames), 'cell') )
    error('lfp_disp:badoptions3', ...
        'Value of ''filenames'' must be a cell array with one string for each filenum' );
end

if ~isempty(evts2avg) && ~(avgflag || ovrflag)
    error('lfp_disp:badoptions4', ...
        '''evtavg'' only works in conjunction with ''avg'' or ''ovr''.');
end

if xcovflag
    for plotnum = 1:numplots
        plotnames{plotnum} = [lfp_FileNames{filenums(2*plotnum - 1)} ...
            ' x ' lfp_FileNames{filenums(2*plotnum)} ];
    end
    axinfo.xlabel = 'Delay, s';
    axinfo.xlim = lfp_Xxcov;
    axinfo.ylim = lfp_Yxcov;
else
    axinfo.xlabel = 'Time, seconds';
    if isempty(window)
        axinfo.xlim = lfp_XLimAll;
    else
        axinfo.xlim = window;
    end
    if isempty(lfp_YLimNum)
        axinfo.ylim = lfp_YLimAll;
    else
        % There is at least one value in lfp_YLimNum, must specify each
        % channel separately in axinfo.ylim
        if isempty(lfp_YLimAll)
            axinfo.ylim = zeros(length(lfp_FileNames),2);
        else
            axinfo.ylim = repmat(lfp_YLimAll,length(lfp_FileNames),1);
        end
        for idx = 1:length(filenums)
            if filenums(idx) <= length(lfp_YLimNum) ...
                    && ~isempty(lfp_YLimNum{filenums(idx)})
                axinfo.ylim(filenums(idx),:) = lfp_YLimNum{filenums(idx)};
            end
        end
    end
end

startedplotting = false; % for Joey's & Theresa's 'eye' plotting to work

% Session names for fig title:
sessionstr = lfp_SessionNames{1};
for k = 2:length(lfp_SessionNames)
    sessionstr = [ sessionstr ' ' lfp_SessionNames{k} ];
end

% Choose proper lfp_multichannel_plot option:
if ovrflag & ~avgflag
    option = 'ovr';
elseif ovrflag & avgflag
    option = 'ovrchan';
elseif avgflag
    option = 'errorcurves';
else
    option = '';
end
if ovrflag
    axinfo.trialnums = [];
end

axinfo.plotflag = plotflag;

%
% Done with argument preprocessing.  Actual analysis/display starts here.
%

% trial loop displays individual trials or collects trial data for
% aggregate processing
for trial = trials
    eventrange = lfp_TrialIndex(trial,1) : lfp_TrialIndex(trial,2);
    
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
            error('lfp_disp:evtbounds1', ...
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
            error('lfp_disp:evtbounds2', ...
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
            find(ismember(trialevents(:,2), lfp_AlignmentRef)), ...
            1 );
    end
    if length(reftime) == 0
        reftime = 0;
        refpoint = 0;
        if norefOKflag
            continue  % abandon whole process of adding row to <trialinfo>
        else
            if marknanbadflag
                % For some reason this condition causes an infinite loop
                % ... Might have fixed that 10-May-2010, did not test -DG
                error('lfp_disp:noref2', ...
                    'Trial %d has no reference event', ...
                    trial );
            else
                warning('lfp_disp:noref2', ...
                    'Trial %d has no reference event', ...
                    trial );
            end
        end
    else
        refisOK(trial) = true;
        if ~multitrigflag
            reftime = reftime(1);
        end
        refpoint = lfp_time2index(reftime);
    end
    
    % <timeinterval>, the time range of data we want to process, is
    % <window> or lfp_XLimAll, limited to the bounds of the trial's
    % recorded time segment.
    if  isempty(evtbounds) && ( (~multitrigflag && reftime == 0) ...
            || (isempty(lfp_XLimAll) && isempty(window)) )
        % We cannot apply <window> or lfp_XLimAll, so we use the whole
        % trial:
        xlimabspoints = [lfp_TrialIndex(trial,3) lfp_TrialIndex(trial,4)];
        timeinterval = [];
    else
        if isempty(window)
            % For single trig, use evtbounds if we have them, or
            % lfp_XLimAll if not:
            if ~multitrigflag && ~isempty(evtbounds)
                timeinterval = ...
                    [lfp_Events(startevtidx,1) lfp_Events(endevtidx,1)] ...
                    - reftime;
                axinfo.xlim = timeinterval;
            else
                if isempty(lfp_XLimAll)
                    error('You should really specify a time window or an lfp_XLimAll');
                else
                    timeinterval = lfp_XLimAll;
                end
            end
        else
            timeinterval = window;
        end
        xlimabs = repmat(timeinterval,size(reftime,1),1) + ...
            repmat(reftime,1,2);
        xlimabspoints = ...
            [lfp_time2index(xlimabs(:,1)) lfp_time2index(xlimabs(:,2))];
    end
    startsample = max(...
        xlimabspoints(:,1), ...
        repmat(lfp_TrialRec(trial,1),size(xlimabspoints,1),1) );
    endsample = min(...
        xlimabspoints(:,2), ...
        repmat(lfp_TrialRec(trial,2),size(xlimabspoints,1),1) );
    if (avgflag || ovrflag)
        % Collect trial info to be used below for computing averaged
        % waveform. <trialinfo> contains the sample indices of
        % the start of selected time interval, end of selected time
        % interval, and reference event.  In case of 'multitrig', multiple
        % rows are added per trial.
        trialinfo = [ trialinfo
            [ startsample endsample refpoint ] ];
        if ovrflag
            axinfo.trialnums = [ axinfo.trialnums
                repmat(trial, size(refpoint)) ];
        end
    elseif eyeflag
        switch lfp_SetupType
            case 'monkey'
                % do the following only once
                if ~startedplotting
                    startedplotting = true;
                    % first flip the colormap 'jet' and save it to a new
                    % variable (for some reason this creates a new figure)
                    mycolormap = colormap(jet);
                    flippedcmap = flipdim(mycolormap,1);
                    hF = gcf;
                    % the following line is prob. redundant
                    set(hF, 'Colormap', flippedcmap);
                    plotcount = 0; % initialize # of plots
                    plothandle = 1; % initialize # of figures
                    a = (find(lfp_SelectedTrials)); % for saving to appropriate filename
                end
                
                if multiwindowflag
                    for trialplot = 1:4
                        plotcount = plotcount + 1; %dumb counter, should change
                        % to subplot exactly the number of trials.
                        % for now, plot up to 16 trials per figure
                        if plotcount > 16
                            if saveflag
                                figfilename = [int2str(lfp_TrialParams{a(1)}(8)) ...
                                    int2str(lfp_TrialParams{a(1)}(9)) ...
                                    int2str(lfp_TrialParams{a(1)}(10)) ...
                                    '_' int2str(plothandle)];
                                saveas(hF, figfilename, 'fig');
                            end
                            
                            if printflag
                                orient(hF, 'landscape');
                                print(hF, '-dwinc');
                            end
                            
                            if saveflag || printflag
                                close(hF);
                            end
                            
                            % then make a new fig
                            hF = figure('Visible', display_str);
                            set(hF, 'Colormap', flippedcmap);
                            plothandle = plothandle + 1;
                            plotcount = 1;
                        end
                        % it would be nice to ensure that the following is
                        % being done in figure(hF), but if we actually use the
                        % call to figure(hF) it will make the figure visible
                        subplot(4,4,plotcount);
                        
                        lfp_plot_eye_jmulti(trialevents, trial, startsample,...
                            endsample, filenums, calflag, trialplot, flippedcmap);
                    end
                else trialplot = 1;
                    plotcount = plotcount + 1; %dumb counter, should change
                    % to subplot exactly the number of trials.
                    % for now, plot up to 16 trials per figure
                    if plotcount > 16
                        if saveflag
                            figfilename = [int2str(lfp_TrialParams{a(1)}(8)) ...
                                int2str(lfp_TrialParams{a(1)}(9)) ...
                                int2str(lfp_TrialParams{a(1)}(10)) ...
                                '_' int2str(plothandle)];
                            saveas(hF, figfilename, 'fig');
                        end
                        
                        if printflag
                            orient(hF, 'landscape');
                            print(hF, '-dwinc');
                        end
                        
                        if saveflag || printflag
                            close(hF);
                        end
                        
                        % then make a new fig
                        hF = figure('Visible', display_str);
                        set(hF, 'Colormap', flippedcmap);
                        plothandle = plothandle + 1;
                        plotcount = 1;
                    end
                    
                    % it would be nice to ensure that the following is
                    % being done in figure(hF), but if we actually use the
                    % call to figure(hF) it will make the figure visible
                    subplot(4,4,plotcount);
                    
                    lfp_plot_eye_jmulti(trialevents, trial, startsample,...
                        endsample, filenums, calflag, trialplot, flippedcmap);
                end
                
                % save or print the last figure
                if (trial == trials(end))
                    if saveflag
                        figfilename = [int2str(lfp_TrialParams{a(1)}(8)) ...
                            int2str(lfp_TrialParams{a(1)}(9)) ...
                            int2str(lfp_TrialParams{a(1)}(10)) ...
                            '_' int2str(plothandle)];
                        saveas(hF, figfilename, 'fig');
                    end
                    if printflag
                        orient(hF, 'landscape');
                        print(hF, '-dwinc');
                    end
                    
                    if saveflag || printflag
                        close(hF);
                    end
                end
            case 'rodent'
                if length(trials) == 1 && ~trackeronlyflag
                    hF = figure('Visible', display_str);
                    data = lfp_plot_eye_rodent(filenums, ...
                        trials, startsample, endsample, argstoeye{:});
                else
                    multitracestartsamples(end+1,1) = startsample;
                    multitraceendsamples(end+1,1) = endsample;
                end
            case 'theresa'
                evts2plot = find( ...
                    lfp_Events(:,1) >= lfp_index2time(startsample) ...
                    & lfp_Events(:,1) <= lfp_index2time(endsample) );
                eyeplotevts = unique( [
                    trialevents
                    lfp_Events(evts2plot,:)
                    ], 'rows');
                % do this only once and only if it is not overlay, in which
                % case it only creates a new figure once
                if ~startedplotting || (~multiflag && ~passoverlayflag)
                    startedplotting = true;
                    % create a figure
                    hF = figure('Visible', display_str);
                    % can use this for mulitple overlay calls
                    %                     hold on
                    %                     hF = figure(100);
                    if (saveflag || printflag) && ~dispflag
                        set(hF, 'Visible', 'off');
                    end
                    plotcount = 0; % initialize # of plots
                    plothandle = 1; % initialize # of figures
                    a = (find(lfp_SelectedTrials)); % for saving to appropriate filename
                    % also for file name use
                    [pathstr, name, ext, versn] = fileparts(lfp_DataDir);
                    % set up to feed multiple arguments to plot_eye_t
                    if passoverlayflag
                        newarg = cell(1,length(argstoeye)+1);
                        % overlay_idx set at very top when first
                        % process 'pass'
                        newarg(1:overlay_idx) = argstoeye(1:overlay_idx);
                        % if the following isn't a string then ismember
                        % complains the next time through
                        newarg{overlay_idx+1} = int2str(length(trials));
                        if length(argstoeye) >= (overlay_idx+1)
                            newarg(overlay_idx+2:end) = argstoeye(overlay_idx+1:end);
                        end
                        argstoeye = newarg;
                    end % end passoverlayflag
                end
                
                if multiflag %w/in 'theresa'
                    plotcount = plotcount + 1; %dumb counter, should change
                    % to subplot exactly the number of trials
                    % for now, plot up to 16 trials per figure
                    if plotcount > (m_plots_t * n_plots_t)
                        % set in lfp_plot_eye_t, not sure needed here
                        % set(hF,'Units','normalized','Position',[0 0 1 1]);
                        if saveflag
                            figfilename = ['Muliti_' pathstr(end-6:end)...
                                '-' name '_' int2str(plothandle) '.fig'];
                            saveas(hF, figfilename, 'fig');
                        end
                        
                        if printflag
                            print(hF, '-dwinc');
                        end
                        
                        if (saveflag || printflag) && ~dispflag
                            close(hF);
                        end
                        
                        % then make a new fig
                        hF = figure('Visible', display_str);
                        if (saveflag || printflag) && ~dispflag
                            set(hF, 'Visible', 'off');
                        end
                        plothandle = plothandle + 1;
                        plotcount = 1;
                    end
                    
                    % it would be nice to ensure that the following is
                    % being done in figure(hF), but if we actually use the
                    % call to figure(hF) it will make the figure visible
                    subplot(m_plots_t, n_plots_t, plotcount)
                    
                    if passflag
                        argstoeye{length(argstoeye)+1} = 'multi';
                        data = lfp_plot_eye_t(eyeplotevts, trial, startsample, ...
                            endsample, filenums, calflag, argstoeye);
                    else
                        data = lfp_plot_eye_t(eyeplotevts, trial, startsample, ...
                            endsample, filenums, calflag, 'multi');
                    end
                    %                 end  %save for plot mult of one trial
                else
                    if passflag
                        data = lfp_plot_eye_t(eyeplotevts, trial, startsample, ...
                            endsample, filenums, calflag, argstoeye);
                    else
                        data = lfp_plot_eye_t(eyeplotevts, trial, startsample, ...
                            endsample, filenums, calflag);
                    end
                end
                
                % save or print the last figure
                if (trial == trials(end)) || ~multiflag
                    if multiflag
                        % set in lfp_plot_eye_t, not sure needed here
                        %set(hF,'Units','normalized','Position',[0 0 1 1]);
                    end
                    if saveflag
                        if ~multiflag
                            uniqueID = lfp_getTrialID(trial);
                            figfilename = ['Trial_' pathstr(end-6:end)...
                                '-' name '-' uniqueID(9:end) '.fig'];
                        else
                            figfilename = ['Muliti_' pathstr(end-6:end)...
                                '-' name '_' int2str(plothandle) '.fig'];
                        end
                        saveas(hF, figfilename, 'fig');
                    end
                    if printflag
                        % orient(hF, 'landscape');
                        print(hF, '-dwinc');
                    end
                    
                    if (saveflag || printflag) && ~dispflag
                        close(hF);
                    end
                end
            case 'simple'
                if isempty(hF)
                    hF = figure;
                    hA = axes('Parent', hF);
                    set(hA, 'NextPlot', 'add');
                    axis(hA, 'equal');
                    set(hA, 'XLim', [0 750]);
                    set(hA, 'YLim', [0 450]);
                    set(hA, 'XDir','normal', 'YDir', 'reverse');
                    lfp_createFigTitle(hA, 'Tracker Plot', trials, [], ...
                        '', '');
                end
                startsamp = lfp_TrialIndex(trial, 3);
                endsamp = lfp_TrialIndex(trial, 4);
                plot(hA, ...
                    lfp_Samples{filenums(1)}(startsamp:33:endsamp), ...
                    lfp_Samples{filenums(2)}(startsamp:33:endsamp), ...
                    'k.', 'MarkerSize', 2);
            otherwise
                error('unrecognized lfp_SetupType');
                
        end
    elseif joyflag  % at this point, ~(avgflag||ovrflag) && ~eyeflag
        switch lfp_SetupType
            case 'monkey'
                lfp_plot_joy(trialevents, trial, startsample, endsample);
            case 'theresa'
                % do nothing
            otherwise
                error('unrecognized lfp_SetupType');
        end
    else
        % => not (avgflag || ovrflag || eyeflag || joyflag), so just plot
        % the trial
        if figflag
            plothandle = trial+figbase;
        elseif multisingflag %TMD added 2/18/08
            if m_plot_num > (m_plots * n_plots)
                m_plot_num = 1;
            end
            if m_plot_num == 1
                figure('Visible', display_str);
            end
            plothandle = subplot(m_plots, n_plots, m_plot_num);
        else
            plothandle = 0;
        end
        if length(startsample) > 1 || length(endsample) > 1
            error('lfp_disp:multirange', ...
                'Multiple start and end indices; bad use of multitrig?');
        end
        trialsamplerange = startsample : endsample;
        if xcovflag
            xcovs = [];
            for plotnum = 1:numplots
                if normflag
                    xcovs(:, plotnum) = dg_xcov(  ...
                        lfp_Samples{filenums(2*plotnum-1)}(trialsamplerange), ...
                        lfp_Samples{filenums(2*plotnum)}(trialsamplerange) );
                else
                    xcovs(:, plotnum) = xcov( ...
                        lfp_Samples{filenums(2*plotnum-1)}(trialsamplerange), ...
                        lfp_Samples{filenums(2*plotnum)}(trialsamplerange), ...
                        'unbiased' )';
                end
            end
            titlestr = sprintf('%s align=%s xcov %s trial %s', ...
                sessionstr, mat2str(lfp_AlignmentRef), ...
                mat2str(timeinterval), lfp_getTrialID(trial) );
            reftime2 = (length(trialsamplerange) - 1) * lfp_SamplePeriod;
            [hF, data] = lfp_multichannel_plot(plothandle, titlestr, ...
                plotnames, xcovs, [], ...
                reftime2, axinfo, option);
        else    % => NOT xcovflag
            % STD mode:
            evts2plot = find( ...
                lfp_Events(:,1) >= lfp_index2time(startsample) ...
                & lfp_Events(:,1) <= lfp_index2time(endsample) );
            [hF, data] = lfp_multichannel_plot(plothandle, ...
                sprintf('align=%s trial %s (%d)', ...
                mat2str(lfp_AlignmentRef), ...
                lfp_getTrialID(trial), trial), ...
                filenums, trialsamplerange, evts2plot, reftime, ...
                axinfo, option, bigevts); %TMD added (trial) 4/10/08
            if multisingflag %TMD 2/18/08
                m_plot_num = m_plot_num + 1;
            end
        end
    end
    if ~(avgflag || ovrflag || isequal(lfp_SetupType,'monkey') ||...
            isequal(lfp_SetupType,'theresa')) && printflag
        print;
    end
end
% End of trial loop

if eyeflag && trackeronlyflag
    hF = figure;
    numtraces = size(multitracestartsamples,1);
    smoothlength=5;
    sampleskip=12;
    calibrated_samples1c = [];
    calibrated_samples2c = [];
    for tracenum = 1:numtraces
        calibrated_samples1c = [ calibrated_samples1c
            reshape(smooth( lfp_Samples{filenums(1)}( ...
            multitracestartsamples(tracenum):sampleskip:multitraceendsamples(tracenum)), ...
            smoothlength ), [], 1) ];
        calibrated_samples2c = [ calibrated_samples2c
            reshape(smooth( lfp_Samples{filenums(2)}( ...
            multitracestartsamples(tracenum):sampleskip:multitraceendsamples(tracenum)), ...
            smoothlength ), [], 1) ];
    end
    plot(calibrated_samples1c, calibrated_samples2c, '.', ...
        'Color', 'k', 'MarkerSize', 4);
    if numtraces > 1
        titlestr = [ sprintf('%s align=%s\n ', ...
            sessionstr, mat2str(lfp_AlignmentRef)) ...
            lfp_getTrialsLabel(trials, trialstyle) ...
            sprintf(' n=%d', length(lfp_enabledTrials(trials))) ];
        title(titlestr, 'Interpreter', 'none');
    else
        titlestr = sprintf('%s align=%s trial %s', ...
            sessionstr, mat2str(lfp_AlignmentRef), ...
            lfp_getTrialID(trials) );
        title(titlestr, 'Interpreter', 'none');
    end
elseif eyeflag && isequal(lfp_SetupType, 'rodent') ...
        && ~isempty(multitracestartsamples)
    hF = figure('Visible', display_str);
    data = lfp_plot_eye_rodent(filenums, ...
        trials, multitracestartsamples, multitraceendsamples, argstoeye{:});
end

% if (avgflag || ovrflag), <trialinfo> now contains pointers to the data to
% be processed.

if avgflag || ovrflag
    if isempty(window)
        if ~isempty(lfp_XLimAll)
            xlimspec = [];
        else
            xlimspec = lfp_XLimAll;
        end
    else
        xlimspec = window;
    end
    if notruncflag && ~isempty(xlimspec)
        xlimallpts = round(xlimspec/lfp_SamplePeriod);
        istoonarrow = trialinfo(:,1) - trialinfo(:,3) >  xlimallpts(1) ...
            | trialinfo(:,2) - trialinfo(:,3) < xlimallpts(2);
        if any(istoonarrow)
            warning( 'lfp_disp:toonarrow', ...
                'Dropping %d triggers due to notrunc option', ...
                sum(istoonarrow) );
            trialinfo(istoonarrow, :) = [];
        end
    end
    if ~isempty(trialinfo) && any(trialinfo(:,2) - trialinfo(:,1) < 2)
        error('lfp_disp:nodata2', ...
            'Time interval to display has collapsed to a single point');
    end
    ntrigs = size(trialinfo,1);
    if ntrigs == 0
        if norefOKflag
            warning('lfp_disp:notrigs', ...
                'No events were found for lfp_AlignmentRef = %s', ...
                dg_thing2str(lfp_AlignmentRef));
            return
        else
            error('lfp_disp:notrigs', ...
                'No events were found for lfp_AlignmentRef = %s', ...
                dg_thing2str(lfp_AlignmentRef));
        end
    end
    % Abort if trials without ref event were included
    if any(trialinfo(:,3)==0)
        badtrials = find(trialinfo(:,3)==0);
        error('lfp_disp:noRefEvent', ...
            'Trial(s) %s had no AlignmentRef event', ...
            dg_canonicalSeries(trials(badtrials)));
    end
    % find maximum time range present in all trials
    pointsbefore = min(trialinfo(:,3) - trialinfo(:,1));
    pointsafter = min(trialinfo(:,2) - trialinfo(:,3));
    % compute avgwave (averaged waveform) for each filenum, one column per
    % file
    avgwave = zeros(pointsbefore + pointsafter + 1, length(filenums));
    if avgxcovflag || ovrflag
        samplematrix = zeros( ...
            size(trialinfo, 1), ...
            pointsbefore+pointsafter+1, ...
            length(filenums) );
    end
    evtmatrix = [];
    trialmatrix = []; % for debugging
    for filenumidx = 1:length(filenums)
        channel = filenums(filenumidx);
        if rectflag
            samples = abs(lfp_Samples{channel});
        else samples = lfp_Samples{channel};
        end
        if avgxcovflag || ovrflag
            % save the raw samples for each channel for each trial in
            % samplematrix(trialidx, sampleidx, filenumidx)
            for trialidx = 1:size(trialinfo, 1)
                samplematrix(trialidx, :, filenumidx) = ...
                    samples( trialinfo(trialidx,3) - pointsbefore ...
                    : trialinfo(trialidx,3) + pointsafter );
            end
            if rmEPflag
                EPs = mean(samplematrix, 1);
                samplematrix = samplematrix - repmat(EPs, ...
                    size(samplematrix, 1), 1);
            end
        end
        if ~(avgxcovflag || ovrflag) || ~isempty(evts2avg)
            % compute the average over trials for this channel, in
            % batches if necessary, using samplematrix2 and subtotwave as
            % working storage
            numtrials = size(trialinfo, 1);
            starttrial = 1;
            batchnum = 1;
            numsamples = pointsbefore + pointsafter + 1;
            batchsize = fix(2^20/numsamples);
            numtrialsleft = numtrials;
            trialhasnan = false(0,1);
            while numtrialsleft > 0
                if err4flag
                    thisbatchsize = numtrialsleft;
                else
                    thisbatchsize = min(numtrialsleft, batchsize);
                end
                samplematrix2 = ...
                    zeros(thisbatchsize, numsamples);   % trials x samples
                rownum = 1; % row number in this batch
                % trigidx is index into trialinfo
                for trigidx = starttrial:(starttrial + thisbatchsize - 1)
                    samplematrix2(rownum,:) = ...
                        samples(trialinfo(trigidx,3) - pointsbefore ...
                        : trialinfo(trigidx,3) + pointsafter) ;
                    if ~isempty(evts2avg)
                        reftime = lfp_index2time(trialinfo(trigidx,3));
                        trial = lfp_time2trial(reftime);
                        eventrange = lfp_TrialIndex(trial,1) : ...
                            lfp_TrialIndex(trial,2);
                        if evtavg2flag
                            starttime = ...
                                lfp_index2time(trialinfo(trigidx,3) ...
                                - pointsbefore);
                            endtime = ...
                                lfp_index2time(trialinfo(trigidx,3) ...
                                + pointsafter);
                            eventrange = eventrange( ...
                                lfp_Events(eventrange,1) >= starttime ...
                                & lfp_Events(eventrange,1) <= endtime );
                        end
                        evtsinrange = lfp_Events(eventrange,:);
                        evts2include = ismember(evtsinrange(:,2), evts2avg);
                        evtmatrix = [evtmatrix
                            [evtsinrange(evts2include,1)-reftime ...
                            evtsinrange(evts2include,2)]
                            ];
                        trialmatrix = [trialmatrix; trigidx trial];
                    end
                    rownum = rownum + 1;
                end
                subtotwave(:,batchnum) = sum(samplematrix2,1)';
                if marknanbadflag
                    trialhasnan = [ trialhasnan
                        any(isnan(samplematrix2), 2) ];
                end
                batchnum = batchnum + 1;
                numtrialsleft = numtrialsleft - thisbatchsize;
                starttrial = starttrial + thisbatchsize;
            end
            if err4flag
                % need to loop over all sample points!
                for k = 1:size(samplematrix2,2)
                    ci(1:2,k) = bootci(numboots, @mean, samplematrix2(:,k));
                end
            end
            clear samplematrix2;
            avgwave(:,filenumidx) = sum(subtotwave,2)/numtrials;
            clear subtotwave;
            
            if errflag || err2flag || err4flag
                % compute the std dev over trials for this channel, in
                % batches if necessary
                numtrials = size(trialinfo, 1);
                starttrial = 1;
                batchnum = 1;
                numsamples = pointsbefore + pointsafter + 1;
                batchsize = fix(2^20/numsamples);
                numtrialsleft = numtrials;
                while numtrialsleft > 0
                    thisbatchsize = min(numtrialsleft, batchsize);
                    diffmatrix = ...
                        zeros(thisbatchsize, numsamples);
                    rownum = 1;
                    for trialidx = starttrial:(starttrial + thisbatchsize - 1)
                        diffmatrix(rownum,:) = ...
                            reshape(samples( ...
                            trialinfo(trialidx,3) - pointsbefore ...
                            : trialinfo(trialidx,3) + pointsafter ), ...
                            1, [] ) ...
                            - reshape(avgwave(:,filenumidx), 1, []);
                        rownum = rownum + 1;
                    end
                    subtotdiffsq(:,batchnum) = sum(diffmatrix.^2,1)';
                    batchnum = batchnum + 1;
                    numtrialsleft = numtrialsleft - thisbatchsize;
                    starttrial = starttrial + thisbatchsize;
                end
                clear diffmatrix;
                stdwave = sqrt(sum(subtotdiffsq,2)/(numtrials-1));
                clear subtotdiffsq;
                if errflag
                    avgwave(:,filenumidx,2) = stdwave;
                else % err2flag
                    avgwave(:,filenumidx,2) = ...
                        2*stdwave/sqrt(ntrigs);
                end
            end
            if err3flag
                avgwave(:,filenumidx,2) = dg_findConfLimits(samplematrix)';
            end
            % this assumes SYMMETRIC confidence limits
            if err4flag
                option = 'asymmetric';
                avgwave(:,filenumidx,2:3) = ci';
            end
            if fileflag
                if isempty(filenames)
                    filename = [];
                else
                    filename = filenames{filenumidx};
                end
                if flipflag
                    dg_save_table(num2cell(trials), [], ...
                        samplematrix(:, ...
                        1:filesampling:size(samplematrix,2), ...
                        :)', filename);
                else
                    dg_save_table([], trials, ...
                        samplematrix(:, ...
                        1:filesampling:size(samplematrix,2), ...
                        :), filename);
                end
            end
        end
        if ovrflag
            % Convert samplematrix(trialidx, sampleidx, filenumidx) =>
            % avgwave(sample, file, trial);
            avgwave = permute(samplematrix, [2 3 1]);
            % Do this last in case there is a prior error:
            lfp_ClickedTrials = [];
            % Add trigtimes:
            axinfo.trigtimes = lfp_index2time(trialinfo(:,3));
        end
        if marknanbadflag
            lfp_BadTrials = union(lfp_BadTrials, trials(trialhasnan));
        end
    end
    % plot the avg waves
    if figflag
        plothandle = 999+figbase;
    else
        plothandle = 0;
    end
    
    
    if multitrigflag
        trigcount = ntrigs;
    else
        trigcount = 0;
    end
    if xcovflag
        % compute actual time interval used in average
        timeinterval = ...
            [ -pointsbefore*lfp_SamplePeriod pointsafter*lfp_SamplePeriod];
        titlestr = [ sprintf('%s align=%s xcov %s\navg trials', ...
            sessionstr, mat2str(lfp_AlignmentRef), ...
            mat2str(timeinterval) ), ...
            sprintf( ' %s n=%d multitrig=%d %s', ...
            lfp_getTrialsLabel(trials, trialstyle), ...
            length(lfp_enabledTrials(trials)), ...
            trigcount, erropt ) ];
        % compute xcovs (array of average xcov column vectors, w/ 3rd D if
        % needed for standard deviations)
        if avgxcovflag || ovrflag
            for plotnum = 1:numplots
                % compute xc (xcov row vector) for each trial
                trial_vector = 1:size(trialinfo, 1);
                if shuffleflag
                    trial_vector_perm = trial_vector;
                    while any(trial_vector_perm == trial_vector)
                        trial_vector_perm = randperm(size(trial_vector,2));
                    end
                    trial_vector = trial_vector_perm;
                end
                if isempty(lfp_Xxcov)
                    maxlags = size(samplematrix,2) - 1;
                else
                    maxlags = min(size(samplematrix,2) - 1, ...
                        round(max(abs(lfp_Xxcov))/lfp_SamplePeriod));
                end
                for trialidx = 1:size(trialinfo, 1)
                    if normflag
                        [xc(trialidx, :), lags] = dg_xcov( ...
                            samplematrix(trialidx, :, 2*plotnum-1), ...
                            samplematrix(trial_vector(trialidx), :, 2*plotnum), ...
                            maxlags);
                    else
                        [xc(trialidx, :), lags] = xcov( ...
                            samplematrix(trialidx, :, 2*plotnum-1), ...
                            samplematrix(trial_vector(trialidx), :, 2*plotnum), ...
                            maxlags, 'unbiased' );
                    end
                end
                if avgxcovflag
                    % compute mean & sd of xc:
                    xcovs(:, plotnum, 1) = mean(xc', 2);
                    if errflag
                        xcovs(:,plotnum,2) = std(xc', 0, 2)';
                    end
                    if err2flag
                        xcovs(:,plotnum,2) = 2*std(xc', 0, 2)'/sqrt(ntrigs);
                    end
                    if err3flag
                        xcovs(:,plotnum,2) = std(xc', 0, 2)';
                    end
                    if fileflag
                        if flipflag
                            dg_save_table(num2cell(trials), [], ...
                                xc(:, 1:filesampling:size(samplematrix,2))');
                        else
                            dg_save_table([], trials, ...
                                xc(:, 1:filesampling:size(samplematrix,2)));
                        end
                    end
                end
                if ovrflag
                    % reshape matrix for overlays: xc(trial, points) =>
                    % xcovs(points, plotnum, trial)
                    xcovs(:, plotnum, :) = xc';
                end
            end
        else
            for plotnum = 1:numplots
                if normflag
                    [xcovs(:, plotnum), lags] = dg_xcov( ...
                        avgwave(:, 2*plotnum-1), ...
                        avgwave(:, 2*plotnum) );
                else
                    [xcovs(:, plotnum), lags] = xcov( ...
                        avgwave(:, 2*plotnum-1), ...
                        avgwave(:, 2*plotnum), 'unbiased' );
                end
            end
        end
        reftime2 = maxlags * lfp_SamplePeriod;
        [hF, data] = lfp_multichannel_plot(plothandle, titlestr, ...
            plotnames, xcovs, [], ...
            reftime2, axinfo, option);
    else % ~xcovflag
        if length(lfp_enabledTrials(trials)) == 1
            % Use single-trial format title if there is only one trial
            titlestr = sprintf('%s align=%s %s trial %s', ...
                sessionstr, mat2str(lfp_AlignmentRef), ...
                mat2str(timeinterval), lfp_getTrialID(trial) );
        else
            titlestr = [ sprintf('%s align=%s\navg ', ...
                sessionstr, mat2str(lfp_AlignmentRef)) ...
                lfp_getTrialsLabel(trials, trialstyle) ...
                sprintf(' n=%d', length(lfp_enabledTrials(trials))) ...
                sprintf(' multitrig=%d %s', trigcount, erropt) ];
        end
        reftime2 = pointsbefore * lfp_SamplePeriod;
        if isempty(avgwave) || (isequal(size(avgwave), [1 1]))
            error('lfp_disp:nodata', ...
                'No data were selected for display; check values\nof lfp_BadTrials, lfp_XLimAll, lfp_SelectedTrials.' );
        end
        plotnames = {};
        for k = 1:length(filenums)
            plotnames{k} = sprintf('%s, %s', lfp_FileNames{filenums(k)}, ...
                lfp_SamplesUnits{filenums(k)});
        end
        % Actual call to plot 'avg' and 'ovr':
        if isequal(option, 'ovrchan')
            auxdata = [];
            reftime2 = reftime - lfp_index2time(refpoint) + reftime2;
        else
            auxdata = [reftime2 lfp_AlignmentRef(1)];
        end
        if dBflag
            [hF, data] = lfp_multichannel_plot(plothandle, titlestr, ...
                plotnames, avgwave, ...
                auxdata, reftime2, axinfo, option, 'dB', dBref);
        else
            [hF, data] = lfp_multichannel_plot(plothandle, titlestr, ...
                plotnames, avgwave, ...
                auxdata, reftime2, axinfo, option);
        end
        if isequal(option, 'ovrchan')
            evts2plot = lfp_Events(:,1) >= lfp_index2time(startsample) ...
                & lfp_Events(:,1) <= lfp_index2time(endsample);
            hA = findobj(hF, 'Type', 'axes');
            for anum = 1:length(hA)
                % Do not plot into the legend axes!
                if isempty(get(hA(anum), 'Tag'))
                    lfp_plotEvtMarkers(hA(anum), lfp_Events(evts2plot, :), ...
                        'reftime', reftime);
                end
            end
        end
        if ovrflag
            data = avgwave;
        end
        if ~isempty(evts2avg)
            axeses = findobj(gcf, 'Type', 'axes');
            oldcurrentaxes = get(hF, 'CurrentAxes');
            for hA = reshape(axeses, 1, [])
                set(hF, 'CurrentAxes', hA);
                hold on;
            end
            if size(evtmatrix,1) == 0
                warning('lfp_disp:events', ...
                    'There were no events selected for display.' );
            else
                evtIDs = unique(evtmatrix(:,2))';
                for ID = evtIDs
                    TS = evtmatrix(evtmatrix(:,2)==ID, 1);
                    medianTS = median(TS);
                    numevents = length(TS);
                    prct = prctile(TS, [5, 95, 25, 75]);
                    eventcolor = [];
                    if ID <= length(lfp_EventColors)
                        eventcolor = lfp_EventColors{ID};
                    end
                    if isempty(eventcolor)
                        eventcolor = lfp_EventDefaultColor;
                    end
                    eventname = '';
                    if ID <= length(lfp_EventNames)
                        eventname = lfp_EventNames{ID};
                    end
                    if isempty(eventname)
                        eventname = '';
                    end
                    colorstr = dg_thing2str(eventcolor);
                    colorstr = regexprep(colorstr, '''', '''''');
                    detailstr = sprintf('ID=%d, "%s", color=%s N=%.0f\\nmedian=%6.4g, mean=%6.4g, SD=%6.4g\\nfirst quartile=%6.4g, last quartile=%6.4g\\n5th percentile=%6.4g, 95th percentile=%6.4g', ...
                        ID, eventname, colorstr, numevents, medianTS, mean(TS), std(TS), prct(3), prct(4), prct(1), prct(2));
                    for hA = reshape(axeses, 1, [])
                        set(hF, 'CurrentAxes', hA);
                        hL = plot([ medianTS medianTS ], ...
                            get(hA, 'YLim'), ...
                            'Color', eventcolor );
                        set(hL, ...
                            'ButtonDownFcn', ...
                            ['fprintf(1,''' detailstr '\n'')'] );
                    end
                end
            end
            set(hF, 'CurrentAxes', oldcurrentaxes); % just in case anyone cares
        end
    end
    if printflag
        print;
    end
end

if markrefOKflag
    lfp_SelectedTrials = refisOK;
end
