function [hF, S_mean, timebinctrs, freqs, plotdata, imfs, hht_f, hht_A, ...
    hht_plot, startpt, jobtable] = lfp_hht(trials, filenum, varargin)
%[hF, S_mean, timebinctrs, freqs, plotdata, imfs, hht_f, hht_A, ...
%     hht_plot, startpt, jobtable] = lfp_hht(trials, filenum)
% Performs the Hilbert-Huang Transform on the CSC data in <filenum> taken
% from <trials> determined in the usual fashion by lfp_AlignmentRef and
% lfp_XLimAll.  Averages over all combinations of trials and values of
% minsift, maxp, and margintime. Frequency markings denote bin centers,
% i.e. the average of the upper and lower bounding freqs, regardless of
% whether they are linearly or logarithmically spaced.
%OUTPUTS
% hF - figure handle to spectrogram-style plot
% S_mean - data array for final spectrogram-style plot
% timebinctrs - the centers of the time bins as specified by the 'tstep'
%   option (default is one sample per bin).
% freqs - nominal "center frequencies" for each freq bin, computed by
%   averaging the upper and lower limits.  This means that when using
%   'logf', the "center frequencies" are actually off-center on the log
%   scale.
% plotdata - all the additional values needed to create plots using
%   lfp_plot_hht.
% imfs - a cell vector containing one element for each "job", i.e. each
%   combination of trial and decomposition parameters.  Each cell contains
%   the set of IMFs for that job as returned by dg_emd but with the last
%   row (containing the "residual junk") deleted.
% hht_f - like <imfs>, but contains instantaneous frequency traces in Hz.
% hht_A - like <imfs>, but contains amplitude traces.
% hht_plot - like <imfs>, but contains the logical vector calculated by
%   'plotinfreqlim' option that is true for the IMFs selected for plotting.
% startpt - an integer vector containing the value of <startpt> that is
%   passed to dg_hht_plotfreqs and dg_hht_plotimfs for each "job", which
%   represents the time of first sample point in that job relative to
%   nominal time 0, in samples.  Note that the time interval beginning with
%   <startpt> is the ENTIRE interval analyzed, and thus INCLUDES
%   <margintime> (if any).
% jobtable - a struct with fields 'midx', 'pidx', 'sidx', 'trialidx' that
%   specifies which combination of the parameters <margintime>, <maxp>,
%   <minsift>, and <trial> respectively, by specifying the indices into the
%   arrays containing those values.  The arrays for <margintime> and
%   <minsift> are sorted in descending order and <maxp> is sorted in
%   ascending order before computation starts, so a value of
%   jobtable(k).midx = 3 means that for job <k>, <margintime> had the third
%   largest value of all values specified (see 'margintime' option>.
%   <trials> does not get sorted, but non-enabled trials are removed before
%   computation starts, so jobtable(k).trial = 7 refers to the seventh
%   ENABLED trial in the trial list.
%OPTIONS
% 'autocolor', range - sets color scale to be the highest value in S_mean
%   down to that value minus <range>.
% 'chatty' - but not so much so as 'verbose'; good for monitoring
%   long-running processes.
% 'dB', bgval - convert power to dB (re: 1 lfp_SamplesUnits{filenum}^2)
%   before blurring and trial averaging (but after IMF averaging); the
%   resulting -Inf values in bins that contain zero power are replaced by
%   <bgval>.  If <bgval> is [], then a value 100 dB below the median power
%   computed over all trials is used.
% 'dBpost' - convert power to dB (re: 1 lfp_SamplesUnits{filenum}^2) after
%   blurring and trial averaging; -Inf values are left in place.
% 'dBprepost' - convert power to dB (re: 1 lfp_SamplesUnits{filenum}^2)
%   after blurring but before trial averaging, exactly like creating one
%   gram for each trial using 'dBpost' and then aggregating them with
%   dg_avgGram(...'sem'...); -Inf values are left in place.
% 'dq' - use Direct Quadrature method to calculate instantaneous
%   frequencies from normalized IMFs returned by dg_AMFMdecomp (default
%   method is Normalized Hilbert Transform) 
% 'fblur', fb - blur the final gram by fb Hz in the frequency
%   direction.  In the case of a logarithmic frequency scale, the blur is a
%   fixed number of bins wide equivalent to <fb> Hz in the middle of the
%   frequency scale.  The smoothing operation has unity gain at zero
%   spatial frequency (i.e. when power is constant over all times and
%   frequencies), but note that this may result in very small values if the
%   instantaneous freq data are sparse.  The smoothing kernel is a normal
%   (Gaussian) distribution function whose total width is equal to 2*<fb>
%   plus one point.  If a two-element value is specified for <fb>, the
%   second element is the width parameter or "sigma" of the Gaussian, i.e.
%   the distance from the peak of the Gaussian where its value is 1/sqrt(e)
%   times (about 61% or -2.2 dB) its peak value.  <fb(2)> is also specified
%   in Hz.  The default width is fb(1)/2.
% 'fsmoothing', freq - smooth the frequency traces calculated from the EMD
%   components with a Hanning window 2*n+1 points wide, where n is the
%   number of samples corresponding to one cycle of <freq> Hz.  Default is
%   n = 2 samples.
% 'imfs', imfs - <imfs> is as returned by this function (lfp_hht).
%   Bypasses computation of the IMFs and uses the values supplied.  It is
%   essential when using this feature that all the "job"-related parameters
%   be exactly the same as when <imfs> was originally computed because
%   otherwise, well, you know, something bad might happen.
% 'logf' - bin frequencies on a log scale (default is linear)
% 'margintime', margintime - specifies an amount of time to add to both
%   ends of the analysis time interval and then trim off of the final
%   result to avoid local end effects (note that there are also non-local
%   end effects far from the endpoints, however).  No checking is done to
%   see if this results in record segment boundaries being crossed; it just
%   stupidly adds adjacent samples, which could result in a crash due to
%   negative sample indices if the very first trial is being analyzed from
%   the very beginning. If <margintime> is a vector, then the spectrogram
%   arrays are calculated for each value and averaged in the same manner as
%   different trials.
% 'maxA', maxA - specifies the maximum amplitude that that is displayed as
%   marker size in the instantaneous freq plot produced by 'plotfreqs'
%   option.  Default is to scale each plot to the maximum amplitude over
%   all IMFs.
% 'maxp', maxp - specifies the maximum power that will be permitted
%   at any point in the mean envelope of an IMF, expressed as a fraction of
%   the average power per point in the current residual waveform from which
%   the next IMF is being sought.  Default = .01.  If <maxp> is a vector,
%   then the spectrogram arrays are calculated for each value and averaged
%   in the same manner as different trials.
% 'maxsift', n - similar to S in Huang et al, but provides an arbitrary
%   quitting criterion in case an IMF that meets other criteria cannot be
%   found (default is 5000), i.e. when the total number of sifts done
%   reaches <n>
% 'methods', methods - <methods> is a cell string vector specifying the
%   interpolation methods for dg_emd in element 1 and dg_hht in element 2,
%   ultimately passed through to Matlab function 'interp1'.  Default is
%   {'spline' 'pchip'}. 
% 'minsift', n - same as S in Huang et al pp 2321-2322: specifies the
%   number of sifts to do after finding an IMF with a constant number of
%   extrema and zero crossings (default is 2).  If <n> is a vector, then
%   the spectrogram arrays are calculated for each value and averaged in
%   the same manner as different trials.
% 'noplot' - too complicated to be understood or even described.
% 'numcores', n - use <n> cores for running multiple trials in parallel;
%   default value is 8.  Setting n=1 bypasses the whole file-based
%   high-level parallelization mechanism, which helps for debugging.
% 'numfreqs', n - creates a spectrogram with <n> frequency bins (default is
%   100)
% 'plotamps', imfnums - for each trial, creates a figure with a plot of
%   instantaneous amplitude vs. time for the IMFs that are listed by number
%   in the integer vector <imfnums>.  If 'plotinfreqlim' is also given,
%   then the value in <imfnums> is ignored, but a dummy value (e.g. 0) must
%   still be supplied as a placeholder for <imfnums>.
% 'plotfreqs' - for each trial, creates a figure with a line plot of
%   instantaneous freq vs. time for all IMFs, after smoothing per
%   'fsmoothing', with clickable reporting of IMF#.
% 'plotfreqamps' - for each trial, creates a scatter plot where amplitude
%   is coded as symbol size.  This plot does not have clickable info
%   because under some conditions Matlab cannot correctly discriminate
%   which 'hggroup' was clicked on.
% 'plotimfs' - for each trial, creates a figure with one plot showing all
%   IMFs, with click info.
% 'plotinamplim', amplim - adds an additional constraint to the IMF
%   selection done by 'plotinfreqlim' such that the amplitude of the IMF
%   must be in the interval <amplim> at the same time as the frequency is
%   in <freqlim>.  <amplim> is expressed as a fraction of the maximum
%   amplitude across all IMFs.
% 'plotinamplim2', amplim - same as 'plotinamplim', except that <amplim>
%   refers to the fraction of each IMF's own maximum amplitude.
% 'plotinfreqlim', freqlim - plots only those imfs whose smoothed
%   frequencies at some point enter the interval defined by <freqlim>.
%   Affects only optional plots, not the spectrogram plot or the returned
%   value of <S_mean>.  Overrides the selection of IMFs given with
%   'plotamps'.  <freqlim> is converted to radians per sample, so floating
%   point roundoff error is possible.
% 'plotimfs2' - for each trial, creates one or more figures with one
%   subplot for each IMF, max four subplots per figure.
% 'tblur', tb - like 'fblur', but operates in the time direction; <tb> is
%   in seconds.
% 'sem', level - adds white contour lines to the average figure showing
%   where the standard error of the mean is equal to <level>.  If <level>
%   contains more than one element, then contour lines are drawn at the
%   remaining values in progressively darker shades of gray.  When combined
%   with 'dB', SEM is calculated from log-scaled power.  With 'dBpost', SEM
%   is calculated in dB as 10*log10( (std(S)+mean(S)) / mean(S) ), where S
%   is the linear power.  Note that in the case of 'dBpost',  the
%   legitimacy of calculating confidence limits as +/- 2 * SEM from the
%   mean depends on the assumption that std(S) << mean(S); if the
%   assumption does not hold, then conf limits have to be calculated as
%   10*log10(mean(S)+2*SEM(S)) and 10*log10(mean(S)-2*SEM(S)), which makes
%   the contour plot display irrelevant.
% 'session', sessionname
%   Supplies a default session name to use when specifying Unique Trial IDs
%   instead of internal trial numbers.
% 'single' - returns S_mean, timebinctrs, freqs, imfs, hht_f, hht_A in
%   single-precision.
% 'tstep', dt - specifies the time represented by 1 point in the final
%   display.  Note that since Matlab plots images independently of screen
%   resolution, this does not necessarily correspond to 1 pixel.  The
%   specified value of <dt> is converted to an integral number of samples.
% 'verbose' - invokes the mode of operation in which its name implicitly
%   states that it will perform.
% 'whiten' - multiply power measurements by f^2 to whiten pink signals

%$Rev: 287 $
%$Date: 2012-10-26 20:14:45 -0400 (Fri, 26 Oct 2012) $
%$Author: dgibson $

% I'm not sure which is uglier, explicitly listing every global used here,
% or explicitly passing all the variables into the functions that now have
% nested definitions.  But lfp_declareGlobals and nesting won't coexist.
global lfp_SamplePeriod lfp_SelectedTrials lfp_XLimAll lfp_Samples ...
    lfp_FreqLim lfp_Events lfp_TrialIndex  lfp_FileNames ...
     lfp_CLimAll lfp_ActiveFilenums

% To keep the "static workspace" happy (i.e. static):
lfp_hhtOneTrial_fedge = []; %#ok<NASGU>
lfp_hhtOneTrial_f = [];
lfp_hhtOneTrial_A = [];
lfp_hhtOneTrial_imf = [];
lfp_hhtOneTrial_result = [];

hF = [];
S_mean = [];
S_std = [];
timebinctrs = [];
freqs = [];
plotdata = [];
imfs = [];
hht_f = [];
hht_A = [];
hht_plot = [];
precomputed_imfs = [];
startpt = [];
jobtable = struct( 'midx', [], 'pidx', [], 'sidx', [], 'trialidx', [] );

amplim = [];
autocolor = [];
bgval = [];
bigevts = {};
chattyflag = false;
dbflag = false;
dbpostflag = false;
dbprepostflag = false;
dt = 0;
dqflag = false;
fblur_raw = 0;
freqlim = [];
fsmoothing = [];
logfflag = false;
margintime = 0;
maxA = [];
maxp = .01;
maxsift = [];
methods = {'spline' 'pchip'};
minsift = 2;
numcores = 8;
numfreqs = 100;
plotamps = [];
plotflag = true;
plotfreqsflag = false;
plotfreqampsflag = false;
plotimfsflag = false;
plotinamplim2flag = false;
semlevel = [];
separateflag = false;
session = '';
singleflag = false;
tblur_raw = 0;
verboseflag = false;
whitenflag = false;
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'autocolor'
            argnum = argnum + 1;
            autocolor = varargin{argnum};
            if ~(isnumeric(bgval) && numel(autocolor)==1)
                error('lfp_hht:autocolor', ...
                    'autocolor value must be a numeric scalar');
            end
        case 'chatty'
            chattyflag = true;
        case 'dB'
            dbflag = true;
            argnum = argnum + 1;
            bgval = varargin{argnum};
            if ~(isempty(bgval) || isnumeric(bgval))
                error('lfp_hht:bgval', ...
                    'bgval for ''dB'' option must be numeric or []');
            end
        case 'dBpost'
            dbpostflag = true;
        case 'dBprepost'
            dbprepostflag = true;
        case 'dq'
            dqflag = true;
        case 'fblur'
            argnum = argnum + 1;
            fblur_raw = varargin{argnum};
            if ~isnumeric(fblur_raw) || ~ismember(length(fblur_raw), [1 2])
                error('lfp_hht:badtb', ...
                    '''fblur'' requires a numeric scalar or two-element argument');
            end
        case 'fsmoothing'
            argnum = argnum + 1;
            fsmoothing = varargin{argnum};
            if ~isnumeric(fsmoothing) || numel(fsmoothing) ~= 1
                error('lfp_hht:badfsmoothing', ...
                    '''fsmoothing'' requires a numeric scalar argument');
            end
            fsmoothing = round(1 / (fsmoothing * lfp_SamplePeriod));
        case 'imfs'
            argnum = argnum + 1;
            precomputed_imfs = varargin{argnum};
        case 'margintime'
            argnum = argnum + 1;
            margintime = varargin{argnum};
        case 'maxsift'
            argnum = argnum + 1;
            maxsift = varargin{argnum};
        case 'maxp'
            argnum = argnum + 1;
            maxp = varargin{argnum};
        case 'methods'
            argnum = argnum + 1;
            methods = varargin{argnum};
        case 'minsift'
            argnum = argnum + 1;
            minsift = varargin{argnum};
        case 'noplot'
            plotflag = false;
        case 'numcores'
            argnum = argnum + 1;
            numcores = varargin{argnum};
            if ~isnumeric(numfreqs)
                error('lfp_hht:badnumcores', ...
                    '''numcores'' requires a numeric argument');
            end
        case 'numfreqs'
            argnum = argnum + 1;
            numfreqs = varargin{argnum};
            if ~isnumeric(numfreqs)
                error('lfp_hht:badnumfreqs', ...
                    '''numfreqs'' requires a numeric argument');
            end
        case 'logf'
            logfflag = true;
        case 'plotamps'
            argnum = argnum + 1;
            if argnum > length(varargin)
                error('lfp_hht:plotamps', ...
                    '''plotamps'' requires a list of IMF numbers');
            end
            plotamps = varargin{argnum};
        case 'plotfreqs'
            plotfreqsflag = true;
        case 'plotfreqamps'
            plotfreqampsflag = true;
        case {'plotimfs', 'plotimfs2'}
            plotimfsflag = true;
            if isequal(varargin{argnum}, 'plotimfs2')
                separateflag = true;
            end
        case {'plotinamplim' 'plotinamplim2'}
            if isequal(varargin{argnum}, 'plotinamplim2')
                plotinamplim2flag = true;
            end
            argnum = argnum + 1;
            amplim = varargin{argnum};
        case 'plotinfreqlim'
            argnum = argnum + 1;
            % immediately convert Hz to radians per sample:
            freqlim = lfp_SamplePeriod*2*pi*varargin{argnum};
        case 'sem'
            argnum = argnum + 1;
            semlevel = varargin{argnum};
        case 'session'
            argnum = argnum + 1;
            session = true;
        case 'single'
            singleflag = varargin{argnum};
        case 'tblur'
            argnum = argnum + 1;
            tblur_raw = varargin{argnum};
            if ~isnumeric(tblur_raw) || ~ismember(length(tblur_raw), [1 2])
                error('lfp_hht:badtb', ...
                    '''tblur'' requires a numeric scalar or two-element argument');
            end
        case 'tstep'
            argnum = argnum + 1;
            dt = varargin{argnum};
            if ~isnumeric(dt) || ~isequal(size(dt), [1 1])
                error('lfp_hht:baddt', ...
                    '''tstep'' requires a numeric scalar argument');
            end
        case 'verbose'
            verboseflag = true;
        case 'whiten'
            whitenflag = true;
        otherwise
            error('lfp_hht:badoption', ...
                ['The option "' dg_thing2str(varargin{argnum}) '" is not recognized.'] );
    end
    argnum = argnum + 1;
end

if isempty(fsmoothing)
    fsmoothing = 2;
end

if dbprepostflag + dbpostflag + dbflag > 1
    warning('lfp_hht:dbopts', ...
        '''dB'' and ''dBpost'' are mutually exclusive; using dBpost');
    dbflag = false;
end

dt = max(1, round(dt/lfp_SamplePeriod));

if isempty(trials)
    trials = 1:length(lfp_SelectedTrials);
end
if ischar(trials)
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

% Find sample ranges
if isempty(lfp_XLimAll)
    [interval, rawtrialinfo] = lfp_findCommonTime(trials);
else
    [interval, rawtrialinfo] = lfp_findCommonTime(trials, 'recseg');
end
if any(rawtrialinfo(:,3)==0)
    error('lfp_hht:noref', ...
        'The reference event is missing from trial(s) %s', ...
        dg_canonicalSeries(trials(rawtrialinfo(:,3)==0)) );
end
if ~isempty(lfp_XLimAll)
    xlimpoints = round(lfp_XLimAll/lfp_SamplePeriod);
    % Set the interval for analysis to xlimpoints clipped by <interval>:
    interval(1) = max(xlimpoints(1), interval(1));
    interval(2) = min(xlimpoints(2), interval(2));
end

% Detrend, analyze
if logfflag
    if isempty(lfp_FreqLim)
        freqedges = logspace( floor(log10(lfp_SamplePeriod*2*pi)), ...
            ceil(log10(pi)), ...
            numfreqs+1 );
    else
        if lfp_FreqLim(1) == 0
            error('lfp_hht:logf0', ...
                'lfp_FreqLim(1) == 0, which cannot be displayed on a log scale');
        end
        freqedges = logspace( ...
            floor(log10(lfp_SamplePeriod*2*pi*lfp_FreqLim(1))), ...
            ceil(log10(lfp_SamplePeriod*2*pi*lfp_FreqLim(2))), ...
            numfreqs+1 );
    end
else
    if isempty(lfp_FreqLim)
        freqedges = linspace(0, pi, numfreqs+1);
    else
        freqedges = linspace(lfp_SamplePeriod*2*pi*lfp_FreqLim(1), ...
            lfp_SamplePeriod*2*pi*lfp_FreqLim(2), numfreqs+1);
    end
end

hht_opts = {'fsmoothing', fsmoothing, 'numfreqs', numfreqs, ...
    'method', methods{2}};
if logfflag
    hht_opts{end+1} = 'logf';
end
if dqflag
    hht_opts{end+1} = 'dq'; 
end
if whitenflag
    hht_opts{end+1} = 'whiten'; 
end
jobsstarted = {};
jobsinprogress = {};
jobsfinished = {};

% Each job is a combination of one value from each of margintime, maxp,
% minsift, trials.  Yes, we could just nest the entire guts of the function
% inside this nest of for-loops, but then the indentation would make the
% code illegible.  The job parameters are sorted so that the longest jobs
% will run first, based on the crude approximation that <trials> makes no
% difference, <maxp> makes the biggest difference, and <minsift> and
% <margintime> make progressively smaller differences.
jobtable(length(trials) * length(minsift) * length(maxp) ...
    * length(margintime), 1) = struct( ...
    'midx', [], 'pidx', [], 'sidx', [], 'trialidx', [] );
jobidx = 1;
[~, maxpidx] = sort(maxp);
[~, minsiftidx] = sort(minsift, 'descend');
[~, marginidx] = sort(margintime, 'descend');
for pidx = maxpidx
    for sidx = minsiftidx
        for midx = marginidx
            for trialidx = 1:length(trials)
                jobtable(jobidx) = struct( 'midx', midx, 'pidx', pidx, ...
                    'sidx', sidx, 'trialidx', trialidx );
                jobidx = jobidx + 1;
            end
        end
    end
end

for jobidx = 1:size(jobtable,1)
    emd_opts = {'maxp', maxp(jobtable(jobidx).pidx), 'method', methods{1}};
    trialidx = jobtable(jobidx).trialidx;
    if ~isempty(maxsift)
        emd_opts(end+1:end+2) = {'maxsift', maxsift};
    end
    if chattyflag
        emd_opts{end+1} = 'chatty';   %#ok<AGROW>
    end
    if verboseflag
        emd_opts{end+1} = 'verbose';   %#ok<AGROW>
    end
    if ~isempty(minsift)
        emd_opts(end+1:end+2) = ...
            {'minsift', minsift(jobtable(jobidx).sidx)};
    end
    marginsamples = ...
        round(margintime(jobtable(jobidx).midx)/lfp_SamplePeriod);
    myinterval = interval + [-marginsamples marginsamples];
    if verboseflag
        lfp_log(sprintf( ...
            'Jobs in progress:\n%s', dg_thing2str(jobsinprogress)));
    end
    while numcores > 1 && length(jobsinprogress) >= (numcores - 1)
        % pause for a moment
        perl('dg_sleep.pl', '1');
        
        % Check to see if any jobsinprogress are done - when the
        % process finishes writing its output, it creates a file named
        % "<idstr>.done". Add any such jobs to jobsfinished and remove
        % from jobsinprogress.
        doneidx = [];
        for jx = 1:length(jobsinprogress)
            idstr = jobsinprogress{jx};
            flagfilename = [idstr '.done'];
            if exist(flagfilename, 'file')
                jobsfinished{end+1} = idstr; %#ok<AGROW>
                doneidx(end+1) = jx; %#ok<AGROW>
            end
        end
        jobsinprogress(doneidx) = []; %#ok<AGROW>
    end
    if isempty(precomputed_imfs)
        idxrange = myinterval(1):myinterval(2);
        % <indices> is in trials x samples format:
        indices = idxrange + rawtrialinfo(trialidx,3);
        if isempty(indices)
            error('lfp_hht:nodata', ...
                ['No samples were selected; note that if lfp_XLimAll is\n' ...
                'empty, <window> is clipped to start and end of trial.'])
        end
        if any(indices < 1 | indices > numel(lfp_Samples{lfp_ActiveFilenums(1)}))
            warning('lfp_hht:badindices', ...
                'Analysis intervals exceed recorded data; check margintime.');
            return
        end
        waveform = lfp_Samples{filenum}(indices);
        waveform = detrend(waveform);
    else
        waveform = precomputed_imfs(jobidx);
    end
    if numcores == 1
        % Finish processing the current job within this loop
        [lfp_hhtOneTrial_result, ~, ...
            lfp_hhtOneTrial_f, lfp_hhtOneTrial_A, lfp_hhtOneTrial_imf] = ...
            lfp_hhtOneTrial(...
            {waveform, emd_opts, hht_opts, dt, ...
            plotfreqsflag || plotfreqampsflag || (~isempty(plotamps)), ...
            freqedges} );
        if margintime(jobtable(jobidx).midx) ~= 0
            marginsamples = ...
                round(margintime(jobtable(jobidx).midx)/lfp_SamplePeriod);
            marginbins = round(marginsamples / dt);
            if exist('S', 'var')
                % The formulation in the 'else' clause is subject to
                % rounding errors that result in size mismatch later on
                lfp_hhtOneTrial_result = ...
                    lfp_hhtOneTrial_result( :, ...
                    (marginbins + 1) : (marginbins + size(S,2)) );
            else
                lfp_hhtOneTrial_result = ...
                    lfp_hhtOneTrial_result( :, ...
                    marginbins + 1 : end - marginbins );
            end
        end
        if dbflag
            Slog = 10*log10(lfp_hhtOneTrial_result);
            S(:,:,jobidx) = Slog; %#ok<AGROW>
        else
            S(:,:,jobidx) = lfp_hhtOneTrial_result; %#ok<AGROW>
        end
        % Collect single trial event data
        get_single_trial_evts;
        
        % Choose/plot freqs, amps, and/or imfs
        choose_IMFs;
        if plotfreqsflag || plotfreqampsflag || ...
                plotimfsflag || ~isempty(plotamps)
            do_optional_plots;
        end
        startpt(jobidx) = myinterval(1); %#ok<AGROW>
        imfs{jobidx} = lfp_hhtOneTrial_imf; %#ok<AGROW>
        hht_f{jobidx} = lfp_hhtOneTrial_f / (lfp_SamplePeriod*2*pi); %#ok<AGROW>
        hht_A{jobidx} = lfp_hhtOneTrial_A; %#ok<AGROW>
        hht_plot{jobidx} = plotimf; %#ok<AGROW>
    else
        % Add another job. Identify each trial's files with this process ID
        % and the trial number, <idstr>.
        idstr = sprintf('lfp_hht%dj%d', dg_pid, jobidx);
        jobsinprogress{end+1} = idstr; %#ok<AGROW>
        jobsstarted{end+1} = idstr; %#ok<AGROW>
        % Leave data for a new process in a .mat file marked
        % with <idstr>
        inputfilename = sprintf('%s_input.mat', idstr);
        logfilename = sprintf('%s.log', idstr);
        oldplotfreqsflag = plotfreqsflag; % the horror, the horror
        plotfreqsflag = ...
            plotfreqsflag || plotfreqampsflag || (~isempty(plotamps));
        save(inputfilename, 'waveform', 'emd_opts', 'hht_opts', ...
            'dt', 'plotfreqsflag', 'freqedges');
        plotfreqsflag = oldplotfreqsflag; % end the horror, the horror
        % Run a new matlab in the background:
        cmdstr = ['matlab -automation -nosplash -nojvm -r "cd(''' pwd ...
            '''); path(pathdef); lfp_hhtOneTrial(''' ...
            idstr '''); exit" -logfile ' logfilename ' & ' ];
        if verboseflag
            fprintf(1, '%s %s\n', datestr(now, 0), '===== ps; lsof =====');
            system('ps -ef; lsof');
            fprintf(1, '%s %s\n', datestr(now, 0), '===== top =====');
            system('top -l 1');
        end
        if verboseflag
            lfp_log(sprintf('Running command >%s<', cmdstr));
        end
        result = system(cmdstr);
        if verboseflag
            lfp_log(sprintf('Finished command >%s< result=%d', ...
                cmdstr, result));
        end
        if result
            warning('lfp_hht:result', ...
                'Failed to run job %s, result=%d', idstr, result);
        end
    end
end
if numcores > 1
    % Check to see if any jobsinprogress are done - when the process
    % finishes writing its output, it creates a file named "<idstr>.done".
    % Add any such jobs to jobsfinished and remove from jobsinprogress.
    doneidx = [];
    for jx = 1:length(jobsinprogress)
        idstr = jobsinprogress{jx};
        flagfilename = [idstr '.done'];
        if exist(flagfilename, 'file')
            jobsfinished{end+1} = idstr; %#ok<AGROW>
            doneidx(end+1) = jx; %#ok<AGROW>
        end
    end
    jobsinprogress(doneidx) = [];
    if verboseflag
        lfp_log(sprintf( ...
            'jobsstarted:\n%s\njobsinprogress:\n%s\njobsfinished:%s\n', ...
            dg_thing2str(jobsstarted), ...
            dg_thing2str(jobsinprogress), ...
            dg_thing2str(jobsfinished) ));
    end
    % Wait until all jobs are finished:
    oldjobsfinished = {};
    loopcount = 0;
    maxloop = 200;
    while length(jobsfinished) < length(jobsstarted)
        % do nothing but pause for a moment
        perl('dg_sleep.pl', '1');
        jobsfinished = {};
        for jx = 1:length(jobsstarted)
            idstr = jobsstarted{jx};
            flagfilename = [idstr '.done'];
            if exist(flagfilename, 'file')
                jobsfinished{end+1} = idstr; %#ok<AGROW>
            end
        end
        if verboseflag
            lfp_log(sprintf( ...
                'jobsstarted:\n%s\njobsfinished:%s\n', ...
                dg_thing2str(jobsstarted), ...
                dg_thing2str(jobsfinished) ));
        end
        if isequal(jobsfinished, oldjobsfinished)
            loopcount = loopcount + 1;
        else
            loopcount = 0;
        end
        msg = sprintf('No new jobs have finished for %d iterations', ...
            loopcount);
        if verboseflag
            lfp_log(msg);
            fprintf(1, '%s %s\n', datestr(now, 0), msg);
        end
        if loopcount > maxloop
            if verboseflag
                fprintf(1, '%s %s\n', datestr(now, 0), '===== ps; lsof =====');
                system('ps -ef; lsof');
                fprintf(1, '%s %s\n', datestr(now, 0), '===== top =====');
                system('top -l 1');
            end
            error('lfp_hht:stuck', ...
                '%s', msg);
        end
    end
    % Process and delete the output files
    for jobidx = 1:length(jobsfinished)
        idstr = jobsfinished{jobidx};
        outputfilename = sprintf('%s_output.mat', idstr);
        if ~exist(outputfilename, 'file')
            files = dir('*_output.mat');
            fnames = '';
            for fidx = 1:length(files)
                fnames = sprintf('%s%s\n', fnames, files(fidx).name);
            end
            warnstr = sprintf('%s does not exist.  Files:\n%sjobsfinished=%s', ...
                outputfilename, fnames, dg_thing2str(jobsfinished));
            warning('lfp_hht:botch', '%s', warnstr);
            lfp_log(warnstr);
        else
            load(outputfilename);
            if margintime(jobtable(jobidx).midx) ~= 0
                marginsamples = ...
                    round(margintime(jobtable(jobidx).midx)/lfp_SamplePeriod);
                marginbins = round(marginsamples / dt);
                if exist('S', 'var')
                    % The formulation in the 'else' clause is subject to
                    % rounding errors that result in size mismatch later on
                    lfp_hhtOneTrial_result = ...
                        lfp_hhtOneTrial_result( :, ...
                        (marginbins + 1) : (marginbins + size(S,2)) );
                else
                    lfp_hhtOneTrial_result = ...
                        lfp_hhtOneTrial_result( :, ...
                        marginbins + 1 : end - marginbins );
                end
            end
            if dbflag
                Slog = 10*log10(lfp_hhtOneTrial_result);
                S(:,:,jobidx) = Slog;  %#ok<AGROW>
            else
                S(:,:,jobidx) = lfp_hhtOneTrial_result;  %#ok<AGROW>
            end
            get_single_trial_evts;
            choose_IMFs;
            marginsamples = ...
                round(margintime(jobtable(jobidx).midx)/lfp_SamplePeriod);
            myinterval = interval + [-marginsamples marginsamples];
            if plotfreqsflag || plotfreqampsflag || ...
                    plotimfsflag || ~isempty(plotamps)
                do_optional_plots;
            end
            if exist([idstr '_err.mat'], 'file')
                warnstr = sprintf( ...
                    'There was an error processing job %s; keeping temp files', ...
                    idstr);
                warning('lfp_hht:err', '%s', warnstr);
                lfp_log(warnstr);
            else
                flagfilename = [idstr '.done'];
                delete(flagfilename);
                inputfilename = sprintf('%s_input.mat', idstr);
                delete(inputfilename);
                delete(outputfilename);
                logfilename = sprintf('%s.log', idstr);
                delete(logfilename);
                if verboseflag
                    lfp_log(sprintf( ...
                        'Deleted temp files:\n%s\n%s\n%s', ...
                        flagfilename, inputfilename, outputfilename));
                end
            end
        end
        startpt(jobidx) = myinterval(1); %#ok<AGROW>
        imfs{jobidx} = lfp_hhtOneTrial_imf; %#ok<AGROW>
        hht_f{jobidx} = lfp_hhtOneTrial_f / (lfp_SamplePeriod*2*pi); %#ok<AGROW>
        hht_A{jobidx} = lfp_hhtOneTrial_A; %#ok<AGROW>
        hht_plot{jobidx} = plotimf; %#ok<AGROW>
    end
end

% S now contains one layer for each trial

clear lfp_hhtOneTrial_result;
% dbflag causes problems with averaging because the mean of -Inf with
% anything is still -Inf.  Therefore, we replace -Inf with a value 100 dB
% below the median finite value:
if dbflag
    S_inf = isinf(S);
    if isempty(bgval)
        bgval = prctile(S(~S_inf), 50) - 100;
    end
    S(S_inf) = bgval;
end

% Calculate freqs in Hz, apply lfp_FreqLim etc.
freqs = (freqedges(1:end-1) + freqedges(2:end)) ...
    / (2 * lfp_SamplePeriod*2*pi);
if ~isempty(lfp_FreqLim)
    is_in_freqlim = freqs >= lfp_FreqLim(1) & freqs <= lfp_FreqLim(2);
    freqs = freqs(is_in_freqlim);
    S = S(is_in_freqlim, :, :);
end

% Do the blurring (if any)
if ~isequal(fblur_raw, 0) || ~isequal(tblur_raw, 0)
    fbinwidth = freqs(floor(length(freqs)/2) + 1) - ...
        freqs(floor(length(freqs)/2));
    fblurpts = round(fblur_raw(1)/fbinwidth);
    if length(fblur_raw) == 1
        f_sig = fblur_raw/(2*fbinwidth);
    else
        f_sig = fblur_raw(2)/fbinwidth;
    end
    if isequal(fblur_raw, 0)
        f_win = 1;
    elseif f_sig == 0
        f_win = zeros(2*fblurpts+1, 1);
        f_win(fblurpts+1) = 1; 
    else
        f_win = reshape(normpdf(-fblurpts:fblurpts, 0, f_sig), [], 1);
    end
    tblurpts = round(tblur_raw(1)/(dt*lfp_SamplePeriod));
    if length(tblur_raw) == 1
        t_sig = tblur_raw/(2*dt*lfp_SamplePeriod);
    else
        t_sig = tblur_raw(2)/(dt*lfp_SamplePeriod);
    end
    if isequal(tblur_raw, 0)
        t_win = 1;
    elseif t_sig == 0
        t_win = zeros(2*tblurpts+1, 1);
        t_win(tblurpts+1) = 1; 
    else
        t_win = reshape(normpdf(-tblurpts:tblurpts, 0, t_sig), [], 1);
    end
    hw2 = repmat(f_win, 1, 2*tblurpts+1) .* ...
        repmat(t_win', 2*fblurpts+1, 1);
    % In the case of dbflag, bgval = the minimum value of S, and is not
    % generally close to zero, so it is not acceptable for conv2 to pad
    % with zeros at the edges; therefore we do our own padding.  However,
    % if ~dbflag, bgval is still undefined:
    if ~dbflag
        bgval = 0;
    end
    vstrip = bgval(ones(size(S,1)+2*fblurpts, tblurpts));
    hstrip = bgval(ones(fblurpts, size(S,2)));
    % Blurring is very time-consuming, so we avoid doing it unnecessarily.
    % This also needs to be parallelized.
    for jx = 1:size(S,3)
        S2 = conv2([vstrip [hstrip; S(:,:,jx); hstrip] vstrip], ...
            hw2);
        S(:,:,jx) = S2((2*fblurpts+1):(end-2*fblurpts), ...
            (2*tblurpts+1):(end-2*tblurpts)); %#ok<AGROW>
    end
end

if  dbpostflag
    S_mean = mean(S, 3);
    if ~isempty(semlevel)
        S_std = 10*log10( (std(S, 0, 3) + S_mean) ./ S_mean);
    end
    S_mean = 10*log10(S_mean);
elseif dbprepostflag
    S = 10*log10(S);
    S_mean = mean(S, 3);
    if ~isempty(semlevel)
        S_std = std(S, 0, 3);
    end
else
    % This covers the cases of dbflag and no flag:
    S_mean = mean(S, 3);
    if ~isempty(semlevel)
        S_std = std(S,0,3);
    end
end

tvals = lfp_SamplePeriod * (dt*(0:size(S_mean,2)) + interval(1));
timebinctrs = (tvals(1:end-1) + tvals(2:end)) / 2;

if singleflag
    S_mean = single(S_mean);
    timebinctrs = single(timebinctrs);
    freqs = single(freqs);
    for jobnum = 1:length(imfs)
        imfs{jobnum} = single(imfs{jobnum}); %#ok<AGROW>
        hht_f{jobnum} = single(hht_f{jobnum}); %#ok<AGROW>
        hht_A{jobnum} = single(hht_A{jobnum}); %#ok<AGROW>
    end
end

if isempty(varargin)
    plotdata.optstr = '(defaults only)';
else
    for argix=1:length(varargin)
        if argix == 1
            plotdata.optstr = dg_thing2str(varargin{1});
        else
            plotdata.optstr = [plotdata.optstr ', ' dg_thing2str(varargin{argix})];
        end
    end
end
plotdata.dbflag = dbflag;
plotdata.dbpostflag = dbpostflag;
plotdata.dbprepostflag = dbprepostflag;
plotdata.autocolor = autocolor;
plotdata.trialevts = trialevts;
plotdata.filename = lfp_FileNames{filenum};
plotdata.clim = lfp_CLimAll;
plotdata.trials = trials;
plotdata.xlim = lfp_XLimAll;
plotdata.semlevel = semlevel;
plotdata.evts2mark = evts2mark;
plotdata.reftime = reftime;
plotdata.bigevts = bigevts;
plotdata.S_std = S_std;
plotdata.N = size(S,3);
plotdata.period = lfp_SamplePeriod;
extraline = sprintf('%s Click for option values', ...
    plotdata.filename);
plotdata.figtitle = lfp_createFigTitle([], 'Hilbert-Huang Spectrogram', ...
    plotdata.trials, plotdata.xlim, extraline, plotdata.optstr);


% Plot
if plotflag
    hF = lfp_plot_hht(S_mean, timebinctrs, freqs, plotdata);
end


    function choose_IMFs
        if isempty(freqlim)
            plotimf = true(size(lfp_hhtOneTrial_f,1), 1);
        else
            freqinlim = lfp_hhtOneTrial_f > freqlim(1) ...
                    & lfp_hhtOneTrial_f < freqlim(2);
            if isempty(amplim)
                plotimf = any(freqinlim, 2);
            else
                maxamps = max(lfp_hhtOneTrial_A, [], 2);
                if ~plotinamplim2flag
                    maxamps(:) = max(maxamps);
                end
                ampinlim = lfp_hhtOneTrial_A > repmat( ...
                    amplim(1) * maxamps, 1, size(lfp_hhtOneTrial_A, 2)) ...
                    & lfp_hhtOneTrial_A < repmat( ...
                    amplim(2) * maxamps, 1, size(lfp_hhtOneTrial_A, 2));
                plotimf = any(freqinlim & ampinlim, 2);
            end
        end
    end

    function get_single_trial_evts
        reftime = lfp_index2time(rawtrialinfo(trialidx,3));
        intervaltime = myinterval * lfp_SamplePeriod + reftime;
        trial = trials(trialidx);
        trialevts = lfp_Events( ...
            lfp_TrialIndex(trial,1) : lfp_TrialIndex(trial,2), : );
        evts2mark = trialevts(:,1) >= intervaltime(1) ...
            & trialevts(:,1) <= intervaltime(2);
    end

    function do_optional_plots
        extraline = sprintf('%s Click for option values', ...
            lfp_FileNames{filenum});
        if isempty(varargin)
            plotdata.optstr = '(defaults only)';
        else
            for k=1:length(varargin)
                if k == 1
                    plotdata.optstr = dg_thing2str(varargin{1});
                else
                    plotdata.optstr = [plotdata.optstr ', ' dg_thing2str(varargin{k})];
                end
            end
        end
        if plotfreqsflag || plotfreqampsflag
            hF = dg_hht_plotfreqs(lfp_hhtOneTrial_f, ...
                lfp_hhtOneTrial_A, ...
                lfp_SamplePeriod, myinterval(1), '', ...
                lfp_FreqLim, maxA, [plotfreqsflag plotfreqampsflag], ...
                plotimf);
            for fnum = 1:length(hF)
                hA = findobj(hF(fnum), 'Type', 'axes');
                for anum = 1:length(hA)
                    if ~isequal(get(hA(anum), 'Tag'), 'legend')
                        lfp_plotEvtMarkers( hA(anum), ...
                            trialevts(evts2mark,:), ...
                            'reftime', reftime, ...
                            'bigevts', bigevts );
                        ax2title = hA(anum);
                    end
                end
                lfp_createFigTitle(ax2title, 'Hilbert-Huang Frequencies', trials, ...
                    lfp_XLimAll, extraline, plotdata.optstr);
            end
        end
        if plotimfsflag
            hF = dg_hht_plotimfs(lfp_hhtOneTrial_imf, lfp_SamplePeriod, ...
                myinterval(1), '', separateflag, plotimf);
            for fnum = 1:length(hF)
                hA = findobj(hF(fnum), 'Type', 'axes');
                for anum = 1:length(hA)
                    if ~isequal(get(hA(anum), 'Tag'), 'legend')
                        lfp_plotEvtMarkers( hA(anum), ...
                            trialevts(evts2mark,:), ...
                            'reftime', reftime, ...
                            'bigevts', bigevts );
                        ax2title = hA(anum);
                    end
                end
                if separateflag
                    lfp_createFigTitle(ax2title, 'Hilbert-Huang IMFs', trials, ...
                        lfp_XLimAll, lfp_FileNames{filenum}, plotdata.optstr, 'oneline');
                else
                    lfp_createFigTitle(ax2title, 'Hilbert-Huang IMFs', trials, ...
                        lfp_XLimAll, extraline, plotdata.optstr);
                end
            end
        end
        if ~isempty(plotamps)
            timepts = lfp_SamplePeriod * ( ...
                (1:size(lfp_hhtOneTrial_A,2)) - 1 + myinterval(1) );
            hF = figure;
            hA = axes('Parent', hF, 'NextPlot', 'add');
            colors = get(hA,'ColorOrder');
            legendstr = {};
            if ~isempty(freqlim)
                plotamps = find(plotimf);
            end
            for imfnum = reshape(plotamps, 1, [])
                mycolor = colors(mod(imfnum-1, size(colors,1)) + 1, :);
                hL = plot(hA, timepts, lfp_hhtOneTrial_A(imfnum,:), ...
                    'Color', mycolor);
                set(hL, 'ButtonDownFcn', sprintf('disp(''IMF#%d'');', ...
                    imfnum));
                grid on;
                legendstr{end+1} = sprintf('IMF#%d', imfnum); %#ok<AGROW>
            end
            xlabel('Time, s');
            ylabel('IMF Amplitude');
            legend(hA, legendstr);
            lfp_plotEvtMarkers( hA, ...
                trialevts(evts2mark,:), ...
                'reftime', reftime, ...
                'bigevts', bigevts );
            lfp_createFigTitle(hA, 'Hilbert-Huang Amplitudes', trials, ...
                    lfp_XLimAll, extraline, plotdata.optstr);
        end
    end

end


