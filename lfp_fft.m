function lfp_fft(trials, filenums, window, varargin)
%LFP_FFT shows power spectra and coherence.
%lfp_fft(trials, filenums, window)
%lfp_fft(..., 'avg')
%lfp_fft(..., 'avg', 'err')
%lfp_fft(..., figbase)
%lfp_fft(..., 'file')
%lfp_fft(..., 'fileflip')
%lfp_fft(..., 'ham')
%lfp_fft(..., 'log')
%lfp_fft(..., 'norm')
%lfp_fft(..., 'ovr')
%lfp_fft(..., 'print')
%lfp_fft(..., 'rmdc')
%lfp_fft(..., 'session', sessionname)
%lfp_fft(..., 'showtrialnums')
%lfp_fft(trials, filenums)
%lfp_fft(trials)
%lfp_fft

%lfp_fft(trials, filenums, window)
%  <trials> is as in lfp_disp
%  Plots the FFT computed on time window <window> (which is a 2-element
%  array specifying [start-time end-time]) from the file numbers in
%  <filenums>.  If <filenums> is the empty array or not given, then all
%  files are plotted.  If <window> is empty or not given, it defaults to
%  [0 1]. Time is shown relative to the first
%  lfp_AlignmentRef event between the trial's lfp_NominalTrialStart and
%  lfp_NominalTrialEnd events.  If there is no lfp_AlignmentRef event, or
%  if the trial duration does not include the entire time window, then the
%  trial is skipped with a warning.
%lfp_fft(..., 'avg')
%  Instead of plotting each trial in its own figure window, the data for
%  the trials are averaged after being FFT'd, and then plotted.
%lfp_fft(..., 'avg', 'err')
%  Adds error bars showing the standard deviation computed using the
%  1/(n-1) formula for unbiased estimate of variance (see docs for Matlab's
%  "std" function).
%lfp_fft(..., 'avg', 'file')
%  Adding 'file' to 'avg' causes lfp_fft to save the tabulated per-trial
%  spectra to a spreadsheet at the same point in processing where the 'err'
%  values are calculated.  ('file' can be combined with 'err', too.)  Note
%  that a separate spreadsheet must be saved for each filenum on display.
%lfp_fft(..., 'ovr')
%  Overlays the plots of the specified trials on a single figure.  Note
%  that 'avg' overrides 'ovr' if they are both specified.
%lfp_fft(..., 'print')
%  Same as for lfp_disp.
%lfp_fft(..., 'rmdc')
%  Removes the DC component from wave(s) before computing the coherence or
%  spectrum.  Not exhaustively tested in combination with other options.
%  Unlike some other programs, this removes dc AFTER applying the
%  window, so it's just like setting the power at dc to zero.
%lfp_fft(..., 'ham')
%  Apply a Hamming window to the selected wave data before FFTing.
%lfp_fft(..., 'log')
%  Display in decibels, i.e. 10 * log10 of the power spectrum.
%lfp_fft(..., 'norm')
%  Display the power as proportion of total power; title is marked "N"
%  (default "R" for "raw").  If combined with 'avg', normalization is done
%  after averaging, and if 'err' is also specified, the error bars are
%  normalized with respect to total power rather than total error.  If
%  combined with 'log', it is done before taking log, so 0 dB represents
%  the total signal energy.  If combined with 'log', 'avg', and 'err', then
%  the average and standard deviation are computed first, then normalized,
%  then converted to log units.
%lfp_fft(trials)
%  All files are plotted using the default window.  Note that <trials> can
%  be a single integer to plot just one trial (no square brackets needed).
%lfp_fft
%  Plots all files for all selected trials.  Equivalent to lfp_fft([]).
%lfp_fft(..., figbase)
%  figbase is a number to add to the trial number to yield the figure
%  number.
%lfp_fft(..., 'session', sessionname)
%  Supplies a default session name to use when specifying Unique Trial IDs
%  instead of internal trial numbers.
%lfp_fft(..., 'showtrialnums')
%  Labels figure with trial numbers instead of selection rule

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

lfp_declareGlobals;

if nargin <3 || isempty(window)
    window = [0 1];
elseif ~(strcmp(class(window), 'double') ...
        && isequal(size(window), [1 2]) )
    error('lfp_fft:badwindow', ...
        '<window> must be 1x2 number array.' );
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
figflag = false;
vscale = 'lin';
windowing = 'box';
normflag = false;
errflag = false;
fileflag = false;
flipflag = false;
ovrflag = false;
printflag = false;
rmdcflag = false;
while argnum <= length(varargin)
    if strcmp(class(varargin{argnum}), 'double')
        if figflag
            error('lfp_fft:multifigbase', ...
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
            case 'err'
                errflag = true;
            case 'file'
                avgflag = true;
                trialinfo = [];
                fileflag = true;
            case 'fileflip'
                avgflag = true;
                trialinfo = [];
                fileflag = true;
                flipflag = true;
            case 'lin'
                vscale = 'lin';
            case 'log'
                vscale = 'log';
            case 'ham'
                windowing = 'ham';
            case 'norm'
                normflag = true;
            case 'ovr'
                ovrflag = true;
            case 'print'
                printflag = true;
            case 'rmdc'
                rmdcflag = true;
            case 'session'
                argnum = argnum + 1;
                session = varargin{argnum};
            case 'showtrialnums'
                trialstyle = 'trialnums';
            otherwise
                error('lfp_fft:badoption', ...
                    ['The option "' varargin{argnum} ...
                        '" is not recognized.'] );
        end
    end
    argnum = argnum + 1;
end

if strcmp(class(trials), 'char')
    trials = lfp_parseTrialStr(trials, session);
end
% Apply mask to trials:
trials = lfp_enabledTrials(trials);

if any(~ismember(filenums, lfp_ActiveFilenums))
    error('lfp_fft:noSuchFilenum', ...
        ['You requested non-existent file numbers: ' ...
            num2str(filenums(find( ...
            ~ismember(filenums, lfp_ActiveFilenums) ))) ]);
end

if size(trials,1) > 1
    if size(trials,2) > 1
        error('lfp_fft:badTrials1', '<trials> must be an integer vector.');
    end
    trials = trials';
end
if any(trials > size(lfp_TrialIndex,1) | trials < 1)
    error('lfp_fft:badtrialnum', ...
        [ 'You have specified trials that do not exist: ' ...
            num2str(trials(find( ...
            trials > size(lfp_TrialIndex,1) | trials < 1 ))) ]);
end
if ~(strcmp(class(trials), 'double')) ...
        || ~all(fix(trials) == trials)
    error('lfp_fft:badTrials2', '<trials> must be an integer vector.');
end

% Initialize lots of stuff:
pointsbefore = -round(window(1)/lfp_SamplePeriod);
pointsafter = round(window(2)/lfp_SamplePeriod);
windowlength = pointsbefore + pointsafter + 1;
freqs = (0:floor(windowlength/2))' / (windowlength*lfp_SamplePeriod);
if ~isempty(lfp_FreqLim)
    displayfreqs = find( freqs >= lfp_FreqLim(1) ...
        & freqs <= lfp_FreqLim(2) );
else
    displayfreqs = repmat(true, 1, length(freqs));
end
axinfo.xlabel = 'Frequency, Hz';
axinfo.xlim = [];
trialinfo = [];
% Initialize PowerSpectraArray, which has freqs on dimension 1, channels
% (aka filenums) on dimension 2, and trials on dimension 3.
PowerSpectraArray = ...
    zeros(windowlength, length(filenums), length(trials));
plotnames = lfp_FileNames(filenums);
axinfo.ylim = lfp_YPwrLim;
trialcounter = 0;

for trial = trials
    trialcounter = trialcounter + 1;
    eventrange = lfp_TrialIndex(trial,1) : lfp_TrialIndex(trial,2);

    % find the lfp_AlignmentRef timestamp:
    trialevents = lfp_Events(eventrange,:);
    reftime = trialevents( ...
        find(ismember(trialevents(:,2), lfp_AlignmentRef)), ...
        1 );
    if length(reftime) == 0
        warning('lfp_fft:norefevent', ...
            ['Could not find the lfp_AlignmentRef event in trial ' ...
                num2str(trial) ] );
        continue
    else
        reftime = reftime(1);
    end
    
    trialsamplerange = lfp_TrialIndex(trial,3) : lfp_TrialIndex(trial,4);  
    % Collect trial info to be used below...
    % First find 'refpoint', the absolute sample index that is closest
    % in time to the lfp_AlignmentRef event:
    refpoint = lfp_time2index(reftime);
    % Then add to the trialinfo table, which contains the sample indices of
    % the start of trial, end of trial, and reference event:
    trialinfo = [ trialinfo
        [ lfp_TrialIndex(trial,3) lfp_TrialIndex(trial,4) refpoint ] ];

    if reftime + window(1) < lfp_index2time(lfp_TrialRec(trial,1)) ...
            || reftime + window(2) ...
            > lfp_index2time(lfp_TrialRec(trial,2))
        warning('lfp_fft:windowtoowide', ...
            ['The specified window goes outside the recording of trial ' ...
                num2str(trial) ]);
        continue
    end
    waves = zeros(windowlength, length(filenums));
    samplerange = refpoint - pointsbefore : refpoint + pointsafter;
    switch windowing
        case {'box' 'mt'}
            windowfunc = [];
        case 'ham'
            windowfunc = hamming(length(samplerange));
    end
    for fileidx = 1:length(filenums)
        waves(:, fileidx) = ...
            (lfp_Samples{filenums(fileidx)}(samplerange))';
        if ~isempty(windowfunc)
            waves(:, fileidx) = waves(:, fileidx) .* windowfunc;
        end
    end
    if rmdcflag
        waves = waves - repmat(mean(waves, 1), size(waves, 1), 1);
    end
    spectra = fft(waves);
    PowerSpectraArray(:,:,trialcounter) = spectra .* conj(spectra) ...
        / length(samplerange);
    if ~(avgflag || ovrflag)
        if normflag
            % The factor of 2 is because half the power is in symmetrical
            % components between Nyquist freq and sampling freq.
            PowerSpectraArray = 2 * PowerSpectraArray ./ repmat( ...
                sum(PowerSpectraArray), size(PowerSpectraArray,1), 1);
            normmark = 'N';
        else
            normmark = 'R';
        end
        if strcmp(vscale, 'log')
            PowerSpectraArray = 10 * log10(PowerSpectraArray);
        end
        if figflag
            fignum = figbase + trial;
        else
            fignum = 0;
        end
        lfp_multichannel_plot(fignum, ...
            sprintf( '%s align=%s win=%s %s %s %s trial %s', ...
                lfp_DataDir, mat2str(lfp_AlignmentRef), mat2str(window), ...
                windowing, vscale, normmark, lfp_getTrialID(trial) ), ...
            plotnames, ...
            PowerSpectraArray(displayfreqs, :, trialcounter), ...
            freqs(displayfreqs), 0, axinfo);
        if printflag
            print;
        end
    end
end
if avgflag || ovrflag
    if avgflag
        Spectra2Display = mean(PowerSpectraArray, 3);
        if errflag
            Spectra2Display(:,:,2) = std(PowerSpectraArray, 0, 3);
            option = 'errorcurves';
        else
            option = '';
        end
        if fileflag
            for fileidx = 1:length(filenums)
                if flipflag
                    dg_save_table(num2cell(trials), freqs(displayfreqs), ...
                        permute(squeeze( ...
                        PowerSpectraArray(displayfreqs,fileidx,:) ), ...
                        [2 1] )' );
                else
                    dg_save_table(num2cell(freqs(displayfreqs)), trials, ...
                        permute(squeeze( ...
                        PowerSpectraArray(displayfreqs,fileidx,:) ), ...
                        [2 1] ) );
                end
            end
        end
        aggrmode = 'avg';
    else    % i.e. ovrflag
        Spectra2Display = PowerSpectraArray;
        option = 'ovr';
        aggrmode = 'ovr';
    end
    if normflag
        if avgflag
            % Normalize the avg spectrum and the error spectrum wrt to sum
            % over freqs of the avg spectrum
            Spectra2Display(:,:,:) = 2 * Spectra2Display(:,:,:) ...
                ./ repmat( sum(Spectra2Display(:,:,1), 1), ...
                [size(Spectra2Display,1), 1, size(Spectra2Display,3)]);
        else    % i.e. ovrflag
            % Normalize each spectrum wrt its own sum over frequencies
            Spectra2Display(:,:,:) = 2 * Spectra2Display(:,:,:) ./ repmat( ...
                sum(Spectra2Display(:,:,:), 1), size(Spectra2Display,1), 1);
        end
        normmark = 'N';
    else
        normmark = 'R';
    end
    if strcmp(vscale, 'log')
        Spectra2Display(:,:,:) = 10 * log10(Spectra2Display(:,:,:));
    end
    title = [ sprintf('%s align=%s win=%s %s %s %s\n%s ', ...
            lfp_DataDir, mat2str(lfp_AlignmentRef), mat2str(window), ...
            windowing, vscale, normmark, aggrmode)...
            lfp_getTrialsLabel(trials, trialstyle) ];
    %
    if figflag
        fignum = figbase;
    else
        fignum = 0;
    end
    lfp_multichannel_plot(fignum, ...
        title, ...
        plotnames, ...
        Spectra2Display(displayfreqs, :, :), ...
        freqs(displayfreqs), 0, axinfo, option);
    if printflag
        print;
    end
end