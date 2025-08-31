function hF = lfp_spec(windowing, trials, filenums, moving_win, varargin)
%Shows spectrograms.
%hF = lfp_spec(windowing, trials, filenums, moving_win)
%lfp_spec(..., figbase)
%lfp_spec(..., 'antigray')
%lfp_spec(..., 'arrows')
%lfp_spec(..., 'avg')
%lfp_spec(..., 'bigevts', bigevtcodes)
%lfp_spec(..., 'EP')
%lfp_spec(..., 'evtavg', evts2avg)
%lfp_spec(..., 'hfnorm', f)
%lfp_spec(..., 'k', numtapers)
%lfp_spec(..., 'lin')
%lfp_spec(..., 'log')
%lfp_spec(..., 'logsig')
%lfp_spec('coh', ..., 'minspikes')
%lfp_spec('mt',..., 'multiplot', numrowscols, plotnum)
%lfp_spec(..., 'nocolorbar')
%lfp_spec(..., 'norm')
%lfp_spec(..., 'nw', bandwidth)
%lfp_spec(..., 'p', plevel)
%lfp_spec('coh', ..., 'pop')
%lfp_spec('coh', ..., 'popmag')
%lfp_spec(..., 'print')
%lfp_spec(..., 'printsig')
%lfp_spec(..., 'pad', N)
%lfp_spec(..., 'session', sessionname)
%lfp_spec(..., 'spikenorm')
%lfp_spec(..., 'rmBL', BL)
%lfp_spec(..., 'rmdc')
%lfp_spec(..., 'rmEP')
%lfp_spec(..., 'rotate')
%lfp_spec(..., 'showtrialnums')
%lfp_spec('coh', ..., 'shuffle')
%lfp_spec('coh', ..., 'signifcolor')
%lfp_spec(..., 'norefOK')

%  Returns the handle to the newly created figure.
%  <windowing> is one of:
%   'mt' - multitaper analysis
%   'ham' - Hamming window (single taper analysis)
%   'coh' - multitaper coherogram; in this case <filenums> must contain an
%       even number of elements, and each pair is used to construct one
%       coherogram.
%  <trials> is as in lfp_disp
%  <filenums> is a list of file numbers.  It can be [] to use
%   the default value, which is all selected files.
%  <moving_win> can be a two-element row vector, where the first element is
%   the window width in seconds and the second element is the time change
%   between successive windows in seconds. Alternatively, moving_win can be
%   a single number specifying window width, and the step size is
%   automatically set to one quarter of a window.  Alternatively, it can be
%   [] to use the default.  The default value is 1.  The value on the time
%   axis represents the time at the center of the window; in other words,
%   the time window extends by half the window width on each side of the
%   value on the x axis. If lfp_XLimAll is set, the convention here is to
%   widen the the time range in the data by half the window width
%   (<moving_win(1)>) on each end, so that the final time scale will appear
%   as lfp_XLimAll.
%
%OPTIONS
%  You can specify as many of these as you like.
%lfp_spec(..., figbase)
%  <figbase> is a number to add to (100*trial-number + filenum) to yield
%  the figure number.
%lfp_spec(..., 'antigray')
%  Uses a colormap ranging from white (low values) to black (high values).
%lfp_spec(..., 'arrows')
%  Only works with 'coh'.  Overlay phase arrows on statistically
%  significant regions of coherogram.  WARNING:  when arrows are overlaid,
%  nothing on the underlying plot can be clicked, so event markers will not
%  show details and the magnifying glass will not be usable.
%lfp_spec(..., 'avg')
%  Displays a spectrogram that is the average of the power spectra of the
%  trials specified.  In this case, the only event marker shown is for the
%  lfp_AlignmentRef event.  We make the approximation here that
%  lfp_SamplePeriod is exact and constant, regardless of how many frames of
%  CSC data are included in the time range to be averaged.  If <windowing>
%  is 'coh', then see lfp_coherogram function re: averaging.
%lfp_spec(..., 'bigevts', bigevtcodes)
%  As in lfp_spikeAnalysis.
%lfp_spec(..., 'EP')
%  Averages the waveforms first, then computes spectrum of the average.
%lfp_spec(..., 'evtavg', evts2avg)
%  Only works with 'avg'.  Collects relative event times of events whose
%  IDs are in <evts2avg> and plots an event marker at the median event
%  time.  The clickable info for the event marker includes median, mean,
%  SD, first quartile, last quartile, 5th percentile, 95th percentile.  Due
%  to roundoff errors, results are only good to about 1 part in 1e+4.
%lfp_spec(..., 'hfnorm', f)
%  Works only with 'mt'.  Normalizes each time window to the average power
%  in the band from <f> to half the sample rate.
%lfp_spec(..., 'k', numtapers)
%  Forces number of multitapers to <numtapers>; default is 1.
%lfp_spec(..., 'minspikes')
%  Applies only if <windowing> is 'coh'.  Must be followed by a numeric
%  value. The second of each pair of filenums is considered to be a spike
%  train presented as wave data, i.e. containing no values other than 0 or
%  1.  Checks each time window for each unit, and if any unit has greater
%  than zero but fewer than the specified number of spikes in that time
%  window summed over all trials, then the entire vertical column of the
%  coherogram is colored with the maximum value color.
%lfp_spec('mt',..., 'multiplot', numrowscols, plotnum)
%  Applies only if <windowing> is 'mt' or 'coh', and 'avg' is invoked:
%  creates single subplot in the current figure in the plotnum position.
%  Creates a new figure window ONLY if there are NO figure windows open.
%lfp_spec(..., 'nocolorbar')
%  Suppresses display of colorbar.  Only works for windowing = 'mt' and
%  windowing = 'coh'.
%lfp_spec('coh', ..., 'nw', nw)
%  Forces "time-bandwidth product" setting for multitapers to <nw>; default
%  is 1.8.
%lfp_spec(..., 'nodisplay')
%  To control display visibility. By default the figure will be visible.
%  Specifying 'nodisplay' still creates the figure, but with 'Visible'
%  set to 'off'.
%lfp_spec(..., 'p', plevel)
%  Sets the p level for computing the significance threshold in
%  coherograms.  Default value is .01.
%lfp_spec(..., 'pad', N)
%  For Hamming windowing:
%  Pads the waveform data with zeros to lengthen it to the next power of 2
%  that is higher than N times the length of the raw waveform data.  This
%  has the effect of interpolating more frequency points to yield a
%  smoother spectrogram, and also makes the FFT computation run faster.
%  For 'mt' or 'coh' windowing:
%  Unless moving_win(2) happens already to be 2^k samples long for some k,
%  there is always padding, and N controls how many powers of 2 to skip:
%  N=0 pads to the next greater power of 2, N=1 pads to twice that length,
%  N=2 pads to four times that length, etc.  The default value is N=0.
%lfp_spec('mt', ..., 'lin')
%  Sets color scale to linear; does not affect 'coh' or 'ham'.
%lfp_spec('mt', ..., 'log')
%  Sets color scale to log (dB), which is the default; does not affect
%  'coh' or 'ham'.
%lfp_spec(..., 'logsig')
%  Logs the significance level of the number of p<.01 significant cells in
%  coherogram, computed after applying lfp_FreqLim and lfp_XLimAll.  Also
%  logs some other important params: fns (filenums) win nw k minsp
%  (minspikes)
%lfp_spec('mt', ..., 'norm')
%  Normalize the color scale such that the maximum value in the entire gram
%  is 1.0; does not affect 'coh' or 'ham'.
%lfp_spec(..., 'pop')
%  Only works with <windowing>='coh' as of 11/1/04.
%  Computes population average coherogram over all specified pairs of
%  filenums.  This is done independently from, and after, averaging over
%  trials ('avg' is implicitly specified when you specify 'pop').  When
%  'pop' is specified, it is possible to submit multiple values of <trials>
%  so that there is a different set of trials selected for each pair of
%  filenums; in this case <trials> is a cell array that contains exactly
%  one element for each pair of filenums, and trials{K} contains the value
%  of <trials> that is applied to filenums(2*K-1) and filenums(2*K).  When
%  using 'pop', the figure title does not attempt to show all the channel
%  pairs and trials arrays; instead it shows a unique ID string that is
%  referenced in the lfp_lib.log with further details, plus the string
%  'trials1:' followed by the trials array or rule for the first pair of
%  filenums.  'evtavg' and 'bigevts' do nothing with 'pop'.
%lfp_spec(..., 'popmag')
%  Same as 'pop', except the average is computed on the magnitude of the
%  coherence, not the complex coherence; i.e., the phase information is
%  discarded before averaging, preventing destructive interference.
%lfp_spec(..., 'print')
%  Same as for lfp_disp.
%lfp_spec(..., 'printsig')
%  Like 'print', but only prints if there is a significant number of
%  significant bins, using the same p < 0.1 test as for generating the
%  warnings when using 'shuffle'.  Works only if <windowing> is 'coh'. Logs
%  signif level of cell counts, computed after applying lfp_FreqLim and
%  lfp_XLimAll, and if it prints then it logs an additional message.
%lfp_spec(..., 'rmBL', BL)
%  Subtracts the baseline spectrum <BL> (as returned by lfp_BLspectrum)
%  from each power spectrum.  Works only with windowing = 'mt'.  When color
%  scale is 'log', the spectra are converted to dB first, then subtracted.
%  If BL spans a different frequency range from lfp_FreqLim, the overlap is
%  used.  If BL is sampled on a different frequency grid from the
%  spectrogram, then BL is interpolated to match using interp1.
%lfp_spec('coh', ..., 'rmdc')
%  Removes the DC component from wave(s) in each time window before
%  computing the coherence or spectrum.
%lfp_spec(..., 'rmEP')
%  Removes the average over trials from each trial in wave(s) before
%  computing the coherence or spectrum.  Only works with 'avg'.
%lfp_spec(..., 'rotate')
%  Works only with 'coh'.  Adds 90 degrees to the phase before drawing
%  arrows, so zero phase is straight up instead of straight to the right.
%lfp_spec(..., 'session', sessionname)
%  Supplies a default session name to use when specifying Unique Trial IDs
%  instead of internal trial numbers.
%lfp_spec(..., 'showtrialnums')
%  Labels figure with trial numbers instead of selection rule
%lfp_spec(..., 'shuffle')
%  Applies only if <windowing> is 'coh' and 'avg' is invoked:  in the
%  second channel, the trials are randomly shuffled, subject to the
%  constraint that after the shuffle no trial may appear in its original
%  position (this is the standard control procedure to ensure that the
%  observed average coherences are not flukes, i.e. they should disappear
%  when this option is invoked).  If p < 0.1 for the number of
%  "significant" bins in a shuffled coherogram, then a warning message is
%  issued and a message written to the log file showing the actual p value
%  calculated from the binomial distribution function.
%lfp_spec(..., 'spikenorm')
%  Works only with 'mt'.  Normalizes each time window to the sum of the
%  sample values in the time window; this is intended for use with spike
%  trains where there is a 1 at each spike time and 0 at all other times
%  (and is of dubious value on other types of waveforms, including
%  "EP-removed" spike trains).  The sum is taken before applying 'rmdc', so
%  there is no interaction between the two options.
%lfp_spec('coh', ..., 'signifcolor')
%  Applies only if <windowing> is 'coh' and 'avg' is invoked:  sets color
%  scale minimum value equal to the significance level.  Leaves maximum
%  value as is, so if lfp_CLimAll is set, then its maximum value still
%  applies.
%lfp_spec(..., 'norefOK')
%  Simply skips any trials trials that do not have a reference event.

% This program is set up in multiple levels either to plot a result or to
% save data to be aggregated and displayed at a higher level.  The
% obligatory portion of the code iterates over filenums (or pairs of them),
% and within that over trials.  For averaging, the obligatory iteration
% over trials saves only the start and end times for each trial, so that
% the common time period shared by all trials can be computed; then
% there is a second pass iterating over trials to compute the actual
% spectrogram.  In the case of population averages, the second pass comes
% after completing both nested loops, whereas for trial averages within a
% single (pair of) filenum(s), the second pass is inside the filenums loop.

%$Rev: 299 $
%$Date: 2013-04-24 16:04:30 -0400 (Wed, 24 Apr 2013) $
%$Author: dgibson $

lfp_spec_start_time = tic;
lfp_declareGlobals;

invocation_id = datestr(now, 30);

if nargin < 4 || isempty(moving_win)
    moving_win = 1;
end
if nargin < 3 || isempty(filenums)
    filenums = sort(lfp_ActiveFilenums);
end
% Apply lfp_SelectedFiles as a mask to filenums:
filenums = filenums(find(lfp_SelectedFiles(filenums)));

if nargin < 2 || isempty(trials)
    trials = 1:size(lfp_TrialIndex,1);
end

session = '';
trialstyle = lfp_TrialStyle;
padfactor = 0;
minspikes = 0;
K = 1;
NW = 1.8;
absflag = false;
antigrayflag = false;
arrowflag = false;
avgflag = false;
bigevts = {};
BL = [];
colorbarflag = true;
EPflag = false;
evts2avg = [];
figflag = false;
logsigflag = false;
multiplotflag = false;
norefOKflag = false;
normflag = false;
normmark = 'R';
p = .01;
padflag = false;
popflag = false;
printflag = false;
printsigflag = false;
ratenorm = 0;
rmdcflag = false;
rmEPflag = false;
rotateflag = false;
shuffleflag = false;
signifcolorflag = false;
spikenormflag = false;
vscale = 'log';
display_str = 'on';
argnum = 1;
while argnum <= length(varargin)
    if strcmp(class(varargin{argnum}), 'double')
        if figflag
            error('lfp_spec:multifigbase', ...
                'You have specified figbase more than once.');
        else
            figflag = true;
            figbase = varargin{argnum};
        end
    else
        switch varargin{argnum}
            case 'antigray'
                antigrayflag = true;
            case 'avg'
                avgflag = true;
            case 'arrows'
                arrowflag = true;
            case 'bigevts'
                argnum = argnum + 1;
                bigevts = varargin{argnum};
                try
                    bigevtIDs = cell2mat(bigevts(:,1));
                catch
                    error('lfp_spikeAnalysis:badbigevts', ...
                        'Value for <bigevts> is badly formatted.' );
                end
            case 'EP'
                EPflag = true;
            case 'evtavg'
                argnum = argnum + 1;
                evts2avg = varargin{argnum};
            case 'hfnorm'
                argnum = argnum + 1;
                ratenorm = varargin{argnum};
                normmark = ['HF' num2str(ratenorm)];
            case 'k'
                argnum = argnum + 1;
                K = varargin{argnum};
            case 'lin'
                if isequal(windowing, 'mt')
                    vscale = 'lin';
                end
            case 'log'
                if isequal(windowing, 'mt')
                    vscale = 'log';
                end
            case 'logsig'
                logsigflag = true;
            case 'minspikes'
                argnum = argnum + 1;
                if ~strcmp(class(varargin{argnum}), 'double')
                    error('lfp_spec:badminspikes', ...
                        'The minspikes specification must be a number');
                end
                minspikes = varargin{argnum};
            case 'multiplot'
                multiplotflag = true;
                argnum = argnum + 1;
                numrowscols = varargin{argnum};
                argnum = argnum + 1;
                plotnum = varargin{argnum};
            case 'nocolorbar'
                colorbarflag = false;
             case 'norefOK'
                norefOKflag = true;   
            case 'norm'
                if isequal(windowing, 'mt')
                    normflag = true;
                    normmark = 'N';
                end
            case 'nw'
                argnum = argnum + 1;
                NW = varargin{argnum};
            case 'p'
                argnum = argnum + 1;
                p = varargin{argnum};
            case 'pad'
                padflag = true;
                argnum = argnum + 1;
                if ~strcmp(class(varargin{argnum}), 'double')
                    error('lfp_spec:badpad', ...
                        'The padding factor must be a number');
                end
                padfactor = varargin{argnum};
            case 'print'
                printflag = true;
            case 'printsig'
                printsigflag = true;
            case 'pop'
                avgflag = true;
                popflag = true;
            case 'popmag'
                popflag = true;
                absflag = true;
            case 'rmdc'
                rmdcflag = true;
            case 'rmBL'
                argnum = argnum + 1;
                BL = varargin{argnum};
            case 'rmEP'
                rmEPflag = true;
            case 'rotate'
                rotateflag = true;
            case 'shuffle'
                shuffleflag = true;
            case 'session'
                argnum = argnum + 1;
                session = varargin{argnum};
            case 'showtrialnums'
                trialstyle = 'trialnums';
            case 'signifcolor'
                signifcolorflag = true;
            case 'spikenorm'
                ratenorm = -1;
                normmark = 'SP';
            case 'nodisplay'
                display_str = 'off';    
            otherwise
                error('lfp_spec:badoption', ...
                    ['The option "' varargin{argnum} ...
                        '" is not recognized.'] );
        end
    end
    argnum = argnum + 1;
end

if nargin < 1
    error('lfp_spec:nomode', 'The <windowing> argument is required.');
end
if any(~ismember(filenums, lfp_ActiveFilenums))
    error('lfp_spec:noSuchFilenum', ...
        ['You requested non-existent file numbers: ' ...
            num2str(filenums(find( ...
            ~ismember(filenums, lfp_ActiveFilenums) ))) ]);
end
if size(moving_win,1) > 1
    error('lfp_spec:badWindow1', '<moving_win> must be a row vector');
end
if size(moving_win,2) > 2 || ~(strcmp(class(moving_win), 'double'))
    error('lfp_spec:badWindow2', ...
        '<moving_win> must contain one or two numbers');
end
if any(moving_win <= 0)
    error('lfp_spec:badWindow3', ...
        '<moving_win> must contain positive, non-zero numbers');
end
if isequal(size(moving_win), [1 1])
    moving_win = [ moving_win moving_win/4 ];
end

if padflag
    padlabel = sprintf('pad %d', padfactor);
else
    padlabel = '';
end

if padflag
    winpoints = round(moving_win / lfp_SamplePeriod);
    exp = 0;
    while 2^exp < winpoints(1) * padfactor
        exp = exp + 1;
    end
    padlength = 2^exp - winpoints(1);
else
    padlength = 0;
end

if strcmp(windowing, 'coh') && (mod(length(filenums), 2) ~= 0)
    error('lfp_spec:badnumfiles', ...
        'You must specify an even number of files for coherence.' );
end

if ~strcmp(windowing, 'coh') && signifcolorflag
    error('lfp_spec:badoption2', ...
        '''signifcolor'' only works with ''coh''' );
end

if ~strcmp(windowing, 'mt') && ratenorm > 0
    error('lfp_spec:badoption4', ...
        '''hfnorm'' only works with ''mt''' );
end

if ~strcmp(windowing, 'mt') && ratenorm == -1
    error('lfp_spec:badoption5', ...
        '''spikenorm'' only works with ''mt''' );
end

if popflag && ~ismember(windowing, {'coh'})
    error('lfp_spec:badoption3', ...
        '''pop'' and ''popmag'' are not implemented for ''%s''', ...
        windowing );
end

if EPflag && avgflag || EPflag && rmEPflag
    error('''EP'' cannot work with ''avg'' or ''rmEP''.');
end

% Test value(s) of <trials>, convert to cell array if necessary
if isequal(class(trials), 'cell')
    if popflag && strcmp(windowing, 'coh')
        if numel(trials) ~= numel(filenums)/2
            error('lfp_spec:badtrials1', ...
                'Cell array <trials> must contain on element for each pair of filenums' );
        end
    else
        error('lfp_spec:badtrials2', ...
            'Cell array <trials> requires ''pop'' and windowing=''coh''');
    end
else
    trials = {trials};
end
% At this point, <trials> is always a cell array
for cellidx = 1:numel(trials)
    if strcmp(class(trials{cellidx}), 'char')
        trials{cellidx} = lfp_parseTrialStr(trials{cellidx}, session);
    end
    validate_trials(trials{cellidx}, cellidx);
    % Apply mask to trials:
    trials{cellidx} = lfp_enabledTrials(trials{cellidx});
end
% At this point, the cells of <trials> always contain numeric arrays

if shuffleflag || printsigflag || logsigflag
    lfp_log(sprintf(...
        'lfp_spec %s fns=%s win=%s nw=%g k=%g minsp=%g rmdc=%g %s %s', ...
        windowing, dg_thing2str(filenums), ...
        dg_thing2str(moving_win), NW, K, minspikes, rmdcflag, ...
        vscale, normmark ));
end
    
channels2 = [];
if strcmp(windowing, 'coh')
    channels2 = filenums(2:2:end);
    filenums = filenums(1:2:end);
end

for fileidx = 1:length(filenums)
    fprintf('lfp_spec started new filenum\n');
    toc(lfp_spec_start_time);
    trialinfo{fileidx} = [];
    channel = filenums(fileidx);
    if ~isempty(channels2)
        chan_2 = channels2(fileidx);
    end
    if popflag
        if numel(trials) > 1
            trials_list = trials{fileidx};
        else
            trials_list = trials{1};
        end
        trials_string = invocation_id;
        % Log the trial selections for this iteration's filenums
        lfp_log(sprintf('%s files %d:%s %d:%s %s', ...
            invocation_id, channel, lfp_FileNames{channel}, ...
            chan_2, lfp_FileNames{chan_2}, ...
            lfp_getTrialsLabel(trials_list, trialstyle) ));
    else
        trials_list = trials{1};
        trials_string = sprintf('%s n=%d', ...
            lfp_getTrialsLabel(trials_list, trialstyle), ...
            length(lfp_enabledTrials(trials_list)) );
    end
    if isempty(trials_list)
        warning('lfp_spec:noTrials', ...
            'No trials selected for filenum %d', filenums(fileidx) );
        continue
    end
    for trial = trials_list
        fprintf('lfp_spec started new trial\n');
        toc(lfp_spec_start_time);
        % find the lfp_AlignmentRef timestamp:
        eventrange = lfp_TrialIndex(trial,1) : lfp_TrialIndex(trial,2);
        trialevents = lfp_Events(eventrange,:);
        reftime = trialevents( ...
            find(ismember(trialevents(:,2), lfp_AlignmentRef)), ...
            1 );
        if length(reftime) == 0
            if norefOKflag
                continue
            else
                error('lfp_spec:noRefEvent', ...
                [ 'Could not find the reference event in trial ' ...
                    num2str(trial) ]);
            end
        else
            reftime = reftime(1);
        end
        
        % Find 'samplerange', the range of CSC data to analyze, and
        % 'dataStartTime', the time at which 'samplerange' starts relative
        % to the reference event.
        trialsamplerange = lfp_TrialIndex(trial,3) : ...
            lfp_TrialIndex(trial,4);
        if isempty(lfp_XLimAll)
            % Do the whole trial if there is no lfp_XLimAll:
            dataStartTime = lfp_index2time(lfp_TrialIndex(trial,3)) - ...
                reftime;
            samplerange = trialsamplerange;
        else
            % Expand the range of data beyond lfp_XLimAll by half a window
            % on each end, but not going beyond the start and end of
            % recording:
            datalim = [ lfp_XLimAll(1) - moving_win(1)/2 ...
                    lfp_XLimAll(2) + moving_win(1)/2 ];
            timerange = reftime + datalim;
            startsample = max(lfp_time2index(timerange(1)), ...
                lfp_TrialRec(trial,1));
            endsample = min(lfp_time2index(timerange(2)), ...
                lfp_TrialRec(trial,2));
            samplerange = startsample : endsample;
            dataStartTime = lfp_index2time(startsample) - reftime;
        end
        
        samp = lfp_Samples{channel};
        % At this point, <samp> contains raw data.  IMPORTANT NOTE: if any
        % signal processing is done to <samp>, then the same must be done
        % to <data2> when running 'coh'.
        infolabel = [ dg_DoubleBackslash(fullfile( ...
            lfp_DataDir, lfp_FileNames{channel}) ) ...
            ' align=' mat2str(lfp_AlignmentRef) ...
            ' win=' mat2str(moving_win(1:2)) ];
        
        infolabel = sprintf('%s rmdc=%d rmEP=%d %s', ...
            infolabel, rmdcflag, rmEPflag, normmark);
        if strcmp(windowing, 'coh')
            channelname2 = lfp_FileNames{chan_2};
            plottype = 'wave';
        end
        if avgflag || popflag || EPflag
            % Collect trial info to be used below for computing averaged
            % spectra.
            % First find the sample index that is closest in time to the
            % lfp_AlignmentRef event:
            refpoint = lfp_time2index(reftime);
            if (reftime - lfp_index2time(refpoint-1)) ...
                    < (lfp_index2time(refpoint) - reftime)
                refpoint = refpoint - 1;
            end
            % Add to the trialinfo table, which contains the sample indices
            % of the start of trial, end of trial, and reference event:
            trialinfo{fileidx}(end+1,:) = ...
                [ samplerange(1) samplerange(end) refpoint ];
        else
            % just plot the trial
            if figflag
                fignum = figbase + 100*trial + channel;
                try
                    close(fignum);
                end
                hF = figure(fignum, 'Visible', display_str);
            else
                hF = figure('Visible', display_str) ;
            end
            infolabel = sprintf('%s\nTrial %s', ...
                infolabel, lfp_getTrialID(trial) );
            
            label.start = dataStartTime;
            switch windowing
                case 'mt'
                    freqlim = lfp_FreqLim;
                    if ~isempty(BL)
                        if isempty(freqlim)
                            freqlim(1) = max(0, BL.f(1));
                            freqlim(2) = min(1/(2*lfp_SamplePeriod), BL.f(end));
                        elseif ~isequal(BL.f([1 end]), lfp_FreqLim)
                            freqlim(1) = max(lfp_FreqLim(1), BL.f(1));
                            freqlim(2) = min(lfp_FreqLim(2), BL.f(end));
                        end
                    end
                    [S,winmid,f]=lfp_mtspecgram2( samp(samplerange)',...
                        moving_win,...
                        [NW K],padfactor,1/lfp_SamplePeriod,...
                        freqlim,0,0,rmdcflag,ratenorm );
                    if ~isempty(BL)
                        BLidx = (BL.f >= freqlim(1) & BL.f <= freqlim(2));
                        if length(BL.f(BLidx)) ~= length(f)
                            intrpBL = reshape( interp1( BL.f(BLidx), ...
                                BL.sum(BLidx), f ), [], 1) / BL.N ;
                        else
                            intrpBL = BL.sum(BLidx)/BL.N;
                        end
                    end
                    t = winmid * lfp_SamplePeriod + label.start;
                    if normflag
                        S = S/max(max(S));
                    end
                    if isequal(vscale, 'log')
                        S = 10*log10(S);
                        cbarlabel = 'Power, dB';
                        if ~isempty(BL)
                            intrpBL = 10*log10(intrpBL);
                        end
                    else
                        cbarlabel = 'Power (linear scale)';
                    end
                    if ~isempty(BL)
                        S = S - repmat(intrpBL', [size(S,1) 1 size(S,3)]);
                    end
                    if ~colorbarflag
                        cbarlabel = 'none';
                    end
                    % special labelling for subplots:
                    if multiplotflag
                        xaxlabel = lfp_AlignmentRef;
                        yaxlabel = '';
                        titlestr = '';
                    else
                        xaxlabel = 'Time, s';
                        yaxlabel = 'Frequency, Hz';
                        titlestr = sprintf('%s\nMulti-taper nw=%g k=%g pad=%g', ...
                            infolabel, NW, K, padfactor);
                    end
                    [hI, hCB] = dg_showGram(hF, t, f, S', ...
                        sprintf('%s\nMulti-taper nw=%g k=%g pad=%g', ...
                        infolabel, NW, K, padfactor), ...
                        xaxlabel, yaxlabel, ...
                         cbarlabel);
                case 'ham'
                    label.title = [infolabel ' Hamming ' padlabel];
                    lfp_stspecgram(@hamming, samp(samplerange)', ...
                        moving_win, 1/lfp_SamplePeriod, ...
                        label, padlength, rmdcflag );
                case 'coh'
                    % Note that this <data2> is used only within
                    data2 = lfp_Samples{chan_2}(samplerange);
                    datatype2 = 0;
                    label.title = sprintf('%s\nvs. %s %s nw=%g k=%g pad=%g', infolabel, ...
                        plottype, channelname2, NW, K, padfactor);
                    lfp_coherogram3(label, samp(samplerange), data2, ...
                        moving_win, 0, padfactor, colorbarflag, 0, arrowflag, ...
                        minspikes, NW, K, rmdcflag, rotateflag, .01 );
                    cbarlabel = 'Coherence'; % added by JF 080408
                otherwise
                    error('lfp_spec:nosuchmode', ...
                        ['The specified windowing "' windowing '" is not recognized.'] );
            end
            if ~isempty(lfp_CLimAll)
                switch windowing
                    case 'mt'
                        hCB = dg_recolorGram(hCB, lfp_CLimAll, hI);
                    otherwise
                        caxis(lfp_CLimAll);
                        colorbar;
                end
            end
            %Plot event markers
            hA = gca;
            evts2mark = lfp_TrialIndex(trial,1) : lfp_TrialIndex(trial,2);
            lfp_plotEvtMarkers(hA, lfp_Events(evts2mark,:), ...
                'reftime', reftime, ...
                'bigevts', bigevts);
        end
        if antigrayflag
            antigray = repmat((127/128:-1/128:0)', 1, 3);
            colormap(antigray);
        end
        if ~avgflag && ~popflag && printflag
            print;
        end
    end
    fprintf('lfp_spec finished first trials loop\n');
    toc(lfp_spec_start_time);
    
    if avgflag || popflag || EPflag
        if shuffleflag
            avgtype = 'shuffle';
        elseif EPflag
            avgtype = 'EP';
        else
            avgtype = 'avg';
        end
        infolabel = sprintf('%s\n%s %s', ...
            infolabel, avgtype, trials_string );
        % find maximum time range present in all trials
        pointsbefore = min(trialinfo{fileidx}(:,3) - trialinfo{fileidx}(:,1));
        pointsafter = min(trialinfo{fileidx}(:,2) - trialinfo{fileidx}(:,3));
        label.start = -pointsbefore * lfp_SamplePeriod;
        % construct matrix of wave data (and event data if requested);
        % samplematrix is in trials X samples format.
        samplematrix = [];
        evtmatrix = [];
        samples = lfp_Samples{channel};
        for row = 1:size(trialinfo{fileidx}, 1)
            samplematrix = [ samplematrix
                reshape( samples( trialinfo{fileidx}(row,3) - pointsbefore ...
                : trialinfo{fileidx}(row,3) + pointsafter ), 1, [] )];
            if ~isempty(evts2avg)
                reftime = lfp_index2time(trialinfo{fileidx}(row,3));
                starttime = lfp_index2time(trialinfo{fileidx}(row,3) - pointsbefore);
                endtime = lfp_index2time(trialinfo{fileidx}(row,3) + pointsafter);
                trial = lfp_time2trial(reftime);
                eventrange = lfp_TrialIndex(trial,1) : lfp_TrialIndex(trial,2);
                trialevents = lfp_Events(eventrange,:);
                evts2include = ismember(trialevents(:,2), evts2avg) ;
                evtmatrix = [evtmatrix
                    [trialevents(evts2include,1)-reftime trialevents(evts2include,2)]
                    ];
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % IMPORTANT MODIFICATION NOTE: remember that anything that is done
        % to samplematrix must also be done to samplematrix2!
        if rmEPflag || EPflag
            EP = mean(samplematrix, 1);
            if rmEPflag
                samplematrix = samplematrix ...
                    - repmat(EP, size(samplematrix,1), 1);
            elseif EPflag
                samplematrix = EP;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %Plot the spectrogram
        if ~(isequal(windowing, 'coh') && popflag)
            if figflag
                fignum = figbase + channel;
                try
                    close(fignum);
                end
                hF = figure(fignum, 'Visible', display_str);
            else
                if multiplotflag
                    hAsub = subplot(numrowscols(1), ...
                        numrowscols(2), plotnum );
                    hF = get(hAsub, 'Parent');
                else
                    hF = figure('Visible', display_str);
                end
            end
        end
        if any(any(isnan(samplematrix)))
            error('lfp_spec:gotnan', 'NaN in input data');
        end
        switch windowing
            case 'mt'
                if ~isempty(BL)
                    if isempty(lfp_FreqLim)
                        freqlim = [ BL.f(1)  BL.f(end) ];
                    elseif ~isequal(BL.f([1 end]), lfp_FreqLim)
                        freqlim(1) = max(lfp_FreqLim(1), BL.f(1));
                        freqlim(2) = min(lfp_FreqLim(2), BL.f(end));
                    end
                end
                [S,winmid,f]=lfp_mtspecgram2( samplematrix',...
                    moving_win,...
                    [NW K],padfactor,1/lfp_SamplePeriod,...
                    lfp_FreqLim,0,1,rmdcflag,ratenorm );
                if ~isempty(BL)
                    BLidx = (BL.f >= freqlim(1) & BL.f <= freqlim(2));
                    if length(BL.f(BLidx)) ~= length(f)
                        intrpBL = reshape( interp1( BL.f(BLidx), ...
                                BL.sum(BLidx), f ), [], 1) / BL.N ;
                    else
                        intrpBL = BL.sum(BLidx)/BL.N;
                    end
                end
                t = winmid * lfp_SamplePeriod + label.start;
                if normflag
                    S = S/max(max(S));
                end
                if isequal(vscale, 'log')
                    S = 10*log10(S);
                    cbarlabel = 'Power, dB';
                    if ~isempty(BL)
                        intrpBL = 10*log10(intrpBL);
                    end
                else
                    cbarlabel = 'Power (linear scale)';
                end
                if ~isempty(BL)
                    S = S - repmat(intrpBL', [size(S,1) 1 size(S,3)]);
                end
                if ~colorbarflag
                    cbarlabel = 'none';
                end
                % special labelling for subplots:
                    if multiplotflag
                        xaxlabel = lfp_AlignmentRef;
                        yaxlabel = '';
                        titlestr = '';
                        hGram = hAsub;
                    else
                        xaxlabel = 'Time, s';
                        yaxlabel = 'Frequency, Hz';
                        titlestr = sprintf('%s\nMulti-taper nw=%g k=%g pad=%g', ...
                            infolabel, NW, K, padfactor);
                        hGram = hF;
                    end

                [hI, hCB] = dg_showGram(hGram, t, f, S', titlestr, ...
                    xaxlabel, yaxlabel, ...
                    cbarlabel );
            case 'ham'
                label.title = [infolabel ' Hamming ' padlabel];
                lfp_stspecgram(@hamming, samplematrix', moving_win, ...
                    1/lfp_SamplePeriod, ...
                    label, padlength, rmdcflag );
            case 'coh'
                if popflag
                    % We have to do population average, so at this point
                    % do nothing (just build up trialinfo).
                else
                    % 2nd channel is wave data (presumably lfp_spike2wave
                    % output); construct second matrix of sample data
                    samplematrix2 = [];
                    samples2 = lfp_Samples{chan_2};
                    for row = 1:size(trialinfo{fileidx},1)
                        samplematrix2 = [ samplematrix2
                            reshape( ...
                            samples2( trialinfo{fileidx}(row,3) - pointsbefore ...
                            : trialinfo{fileidx}(row,3) + pointsafter ), ...
                            1, [] ) ];
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % IMPORTANT MODIFICATION NOTE: remember that anything that is done
                    % to samplematrix must also be done to samplematrix2!
                    if rmEPflag || EPflag
                        EP = mean(samplematrix2, 1);
                        if rmEPflag
                            samplematrix2 = samplematrix2 ...
                                - repmat(EP, size(samplematrix2,1), 1);
                        elseif EPflag
                            samplematrix2 = EP;
                        end
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % special labelling for subplots:
                    if multiplotflag
                        label.xax = lfp_AlignmentRef;
                        label.yax = '';
                        label.title = '';
                        label.hF = hF;
                    else
                        label.title = sprintf('%s\nvs. %s %s nw=%g k=%g pad=%g', infolabel, ...
                            plottype, channelname2, NW, K, padfactor);
                    end
                    numrows = size(samplematrix2,1);
                    if (numrows > 1) && shuffleflag
                        % Use a list of random numbers to create a random
                        % sorting order, then sort one of the multi-trial data
                        % matrices by rows in that order.  If any trials are
                        % still in their original positions, repeat.
                        randidx = (1:numrows)';
                        while any(randidx == (1:numrows)')
                            randomlist = rand(numrows,1);
                            [b, randidx] = sort(randomlist);
                        end
                        samplematrix2 = samplematrix2(randidx, :);
                    end
                    [thresh, cohdisp] = lfp_coherogram3(label, ...
                        samplematrix, samplematrix2, ...
                        moving_win, 1, padfactor, colorbarflag, 0, arrowflag, ...
                        minspikes, NW, K, rmdcflag, rotateflag, p );
                    % Compute & log "differential summed coherence", i.e.
                    % the sum of the magnitude of coherence in all cells in
                    % the top half of the coherogram minus that in the
                    % bottom half, divided by the number of cells in each
                    % half.  Set all Infs and Nans to zero first.  Note that
                    % at this point, the coherogram is already limited by
                    % lfp_FreqLim and lfp_XLimAll, but we have to apply
                    % thresh explicitly to get what is shown in the
                    % 'signifcolor' graph.6
                    cohdisp2 = max(abs(cohdisp) - thresh, 0);
                    cohdisp2(find(isnan(cohdisp2))) = 0;
                    cohdisp2(find(isinf(cohdisp2))) = 0;
                    halfheight = floor(size(cohdisp,1)/2);
                    topsum = sum(sum(cohdisp2(end-halfheight+1:end, :)));
                    botsum = sum(sum(cohdisp2(1:halfheight, :)));
                    lfp_log(sprintf('differential summed coherence = %d', ...
                        (topsum - botsum) / (halfheight * size(cohdisp,2)) ));
                    if shuffleflag || printsigflag || logsigflag
                        % compute p2, the probability of getting the
                        % observed number or larger of "significant" cells
                        cohdisp2 = cohdisp .* ~isinf(cohdisp);
                        numsignif = sum(sum(abs(cohdisp2) > thresh));
                        p2 = sum(binopdf(...
                            numsignif:numel(cohdisp), numel(cohdisp), p) );
                        lfp_log(sprintf(...
                            'lfp_spec: significance level of cell counts = %4.2d', ...
                            p2 ));
                        if p2 < 0.1
                            if shuffleflag
                                lfp_log(sprintf(...
                                    'Too many significant cells in %s', ...
                                    label.title ));
                                warning('lfp_spec:signifshuffle', ...
                                    'Too many significant cells (p = %4.2d)', p2 );
                            end
                        end
                    end
                end % if popflag ... else            
            otherwise
                error('lfp_spec:nosuchmode', ...
                    ['The specified windowing "' windowing '" is not recognized.'] );
        end
        
        if ~(isequal(windowing, 'coh') && popflag)
            if ~isempty(lfp_CLimAll) || signifcolorflag
                if signifcolorflag
                    colorlim = caxis;
                    if ~isempty(lfp_CLimAll)
                        colorlim(2) = max(thresh*(1+eps), lfp_CLimAll(2));
                    else
                        colorlim(2) = max(thresh*(1+eps), colorlim(2));
                    end
                    colorlim(1) = thresh;
                else
                    % lfp_CLimAll must be non-empty
                    colorlim = lfp_CLimAll;
                end
                caxis(colorlim);
                if colorbarflag
                    colorbar;
                end
            end
            
            %Plot event marker(s)
            hA = gca;
            lfp_plotEvtMarkers(hA, [0 lfp_AlignmentRef(1)]);
            if ~isempty(evts2avg)
                lfp_plotEvtMarkers(hA, evtmatrix, 'stats', ...
                    'bigevts', bigevts);
            end
            if antigrayflag
                antigray = repmat((127/128:-1/128:0)', 1, 3);
                colormap(antigray);
            end
            if printflag
                print;
            end
            if printsigflag && p2 < 0.1
                print;
                lfp_log('lfp_spec: printed');
            end
        end
    end
    fprintf('lfp_spec finished post-processing\n');
    toc(lfp_spec_start_time);
end
fprintf('lfp_spec finished iterating over filenums\n');
toc(lfp_spec_start_time);

% Plot population-averaged coherogram
if isequal(windowing, 'coh') && (popflag)
    if figflag
        fignum = figbase + 100*trial + channel;
        try
            close(fignum);
        end
        hF = figure(fignum, 'Visible', display_str);
    else
        hF = figure('Visible', display_str);
    end
    % Construct cell array of wave data
    % find maximum time range present in all trials
    for fileidx = 1:numel(trialinfo)
        pointsbefore(fileidx) = min(...
            trialinfo{fileidx}(:,3) - trialinfo{fileidx}(:,1) );
        pointsafter(fileidx) = min(...
            trialinfo{fileidx}(:,2) - trialinfo{fileidx}(:,3) );
    end
    pointsbefore = min(pointsbefore);
    pointsafter = min(pointsafter);
    label.start = -pointsbefore * lfp_SamplePeriod;
    for fileidx = 1:numel(trialinfo)
        cohdata1{fileidx} = [];
        cohdata2{fileidx} = [];
        channel = filenums(fileidx);
        chan_2 = channels2(fileidx);
        samples = lfp_Samples{channel};
        samples2 = lfp_Samples{chan_2};
        for row = 1:size(trialinfo{fileidx},1)
            cohdata1{fileidx} = [ cohdata1{fileidx}
                samples(trialinfo{fileidx}(row,3) - pointsbefore ...
                    : trialinfo{fileidx}(row,3) + pointsafter) ];
            cohdata2{fileidx} = [ cohdata2{fileidx}
                samples2(trialinfo{fileidx}(row,3) - pointsbefore ...
                    : trialinfo{fileidx}(row,3) + pointsafter) ];
        end
        numrows = size(cohdata2{fileidx},1);
        if (numrows > 1) && shuffleflag
            % Use a list of random numbers to create a random
            % sorting order, then sort one of the multi-trial data
            % matrices by rows in that order.  If any trials are
            % still in their original positions, repeat.
            randidx = (1:numrows)';
            while any(randidx == (1:numrows)')
                randomlist = rand(numrows,1);
                [b, randidx] = sort(randomlist);
            end
            cohdata2{fileidx} = cohdata2{fileidx}(randidx, :);
        end
    end
    label.title = sprintf('Population %s %s\ntrials1: %s\nalign=%s win=%s',...
        avgtype, invocation_id, lfp_getTrialsLabel(trials{1},trialstyle), ...
        mat2str(lfp_AlignmentRef), ...
        mat2str(moving_win(1:2)));
    [thresh, cohdisp] = lfp_coherogram3(label, cohdata1, cohdata2, ...
        moving_win, 1, padfactor, colorbarflag, 0, arrowflag, ...
        minspikes, NW, K, rmdcflag, rotateflag, p );
    if shuffleflag || printsigflag || logsigflag
        % compute p2, the probability of getting the
        % observed number or larger of "significant" cells
        cohdisp2 = cohdisp .* ~isinf(cohdisp);
        numsignif = sum(sum(abs(cohdisp2) > thresh));
        p2 = sum(binopdf(...
            numsignif:numel(cohdisp), numel(cohdisp), p) );
        lfp_log(sprintf(...
            'lfp_spec: significance level of cell counts = %4.2d', ...
            p2 ));
        if p2 < 0.1
            if shuffleflag
                lfp_log(sprintf(...
                    'Too many significant cells in %s', ...
                    label.title ));
                warning('lfp_spec:signifshuffle', ...
                    'Too many significant cells (p = %4.2d)', p2 );
            end
        end
    end    
    if ~isempty(lfp_CLimAll) || signifcolorflag
        if signifcolorflag
            colorlim = caxis;
            if ~isempty(lfp_CLimAll)
                colorlim(2) = lfp_CLimAll(2);
            else
                colorlim(2) = max(thresh+eps, colorlim(2));
            end
            caxis([thresh colorlim(2)]);
        else
            % lfp_CLimAll must be non-empty
            caxis(lfp_CLimAll);
        end
        colorbar;
    end
    
    %Plot event marker as in lfp_disp
    hA = gca;
    lfp_plotEvtMarkers(hA, [0 lfp_AlignmentRef(1)]);
    if antigrayflag
        antigray = repmat((127/128:-1/128:0)', 1, 3);
        colormap(antigray);
    end
    if printflag
        print;
    end
    if printsigflag && p2 < 0.1
        print;
        lfp_log('lfp_spec: printed');
    end
end


function validate_trials(trials, cellnum)
% Test value of <trials> for validity.  <cellnum> is used only to make
% error messages more helpful; it should be set to the index number of the
% cell being validated when processing an element from a cell array of
% <trials>. If <cellnum> is 0 or missing, then error messages are presented
% without any cell number included, so this is the value to use when
% processing a single value of <trials>.
lfp_declareGlobals;
if nargin < 2
    cellnum = 0;
end
if cellnum
    celltxt = sprintf('Cell #%d\n', cellnum);
else
    celltxt = '';
end
if size(trials,1) > 1
    if size(trials,2) > 1
        error('lfp_spec:badTrials1', ...
            '%s<trials> must be an integer vector.', celltxt);
    end
    trials = trials';
end
if ~(strcmp(class(trials), 'double')) ...
        || ~all(fix(trials) == trials)
    error('lfp_spec:badTrials2', ...
        '%s<trials> must be an integer vector.', celltxt);
end
if any(trials > size(lfp_TrialIndex,1) | trials < 1)
    error('lfp_disp:badtrialnum', ...
        sprintf('%sYou have specified trials that do not exist: %s', ...
        celltxt, num2str(trials(find( ...
        trials > size(lfp_TrialIndex,1) | trials < 1  )))) ...
    );
end
