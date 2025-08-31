function [trigTS, sameTS] = dg_rmStimArtifact(ctrlfilepath, filenames, ...
    offsets, varargin)
% Linearly interpolates stimulus artifacts for fearless filtering.
%INPUTS
% ctrlfilepath: path to a CSC file that contains data used to control the
%   timing of the interpolations.  This could be a recorded CSC file
%   containing large stimulation artifacts, or a control signal that was
%   triggering the stimulator.  It is thresholded at half of its maximum
%   positive-going excursion from the median to generate trigger
%   timestamps. Uses >, <=, so can trigger on waveform that just touches
%   threshold, but will not miss any crossings.  An error is raised if the
%   file does not contain samples that have the statistics expected of a
%   control channel (see local function ctrlfileStats).
% filenames: cell array of names of files to process.  By default, it is
%   assumed that they are in the same directory specified by
%   <ctrlfilepath>.  File extensions must be either '.mat' or '.ncs'.
% offsets: a two-element vector of times with respect to trigger timestamps
%   specifying the start and end times for the linear interpolation.  May
%   be positive or negative, but second element must be greater than first.
%   Specified in seconds if raw recorded timestamps are in microseconds.
%   Second element also determines a "time-out" period during which
%   additional triggers are ignored.
%OUTPUTS
% A '.mat' file is created whose name is composed of the base name of the
%   input file plus the suffix 'rmstim', with the extension '.mat'.  By
%   default, output files are created in the same directory specified by
%   <ctrlfilepath>.  Note that if the input files are in the same
%   directory, and they also have extension '.mat', then they will be
%   overwritten.
% trigTS: a column vector of timestamps of all detected stimulation
%   triggers.  Since this is done using CSC data, the result is only
%   precise to within one sample period.  To be exact, each <trigTS> is the
%   timestamp of the last sample before the control signal exceeded
%   threshold, which by default is halfway between the maximum sample value
%   in the file and the median sample value.
% sameTS: for each file in <filenames>, true if the recorded frame
%   timestamps are the same as the ones in <ctrlfilepath>.
%OPTIONS
% 'dest', dest - path to the directory to which output files will be
%   written.
% 'noVscale' - produces output in raw AD units instead of scaling according
%   to the .ncs file's 'bitvolts' value.
% 'src', src - path to the directory from which input files will be read.
% 'thresh', thresh - fraction of maximum excursion (measured from median
%   value to maximum value) at which to generate trigger timestamp from a
%   rectangular control signal.  Default = 0.5.  (See local function
%   ctrlfileStats for details concerning non-rectangular waveforms.)
% 'threshval', threshval - predetermined fixed threshold at which to
%   generate trigger timestamp from a rectangular control signal.  This
%   also bypasses all calls to ctrlfileStats for the sake of execution
%   speed.  NOTE: thresholding is done on the control file, and no voltage
%   conversion is done on it, even when it is the same file as the signal
%   file.
% 'timeout', timeout - amount of time after a trigger is detected during
%   which any additional triggers are ignored.  Default = <offset(2)>.
% 'verbose' - produces progress report messages in command window.
%NOTES
% Timestamps are internally represented in the same integer units
% (assumed to be microsec) as in the data file.
%   See also dg_convertStimSession.

%$Rev: 253 $
%$Date: 2016-09-20 21:27:40 -0400 (Tue, 20 Sep 2016) $
%$Author: dgibson $

[ctrldir, ~, ext] = fileparts(ctrlfilepath);
if isempty(ctrldir)
    ctrldir = pwd;
end

dest = ctrldir;
offsets = 1e6 * offsets;
src = ctrldir;
thresh = 0.5;
threshvalopt = [];
timeout = offsets(2);
verboseflag = false;
vscale = true;

argnum = 1;
while argnum <= length(varargin)
    if ischar(varargin{argnum})
        switch varargin{argnum}
            case 'dest'
                argnum = argnum + 1;
                dest = varargin{argnum};
            case 'noVscale'
                vscale = false;
            case 'src'
                argnum = argnum + 1;
                src = varargin{argnum};
            case 'thresh'
                argnum = argnum + 1;
                thresh = varargin{argnum};
            case 'threshval'
                argnum = argnum + 1;
                threshvalopt = varargin{argnum};
            case 'timeout'
                argnum = argnum + 1;
                timeout = varargin{argnum};
            case 'verbose'
                verboseflag = true;
            otherwise
                error('dg_rmStimArtifact:badoption', ...
                    'The option %s is not recognized.', ...
                    dg_thing2str(varargin{argnum}));
        end
    else
        error('dg_rmStimArtifact:badoption2', ...
            'The value %s occurs where an option name was expected', ...
            dg_thing2str(varargin{argnum}));
    end
    argnum = argnum + 1;
end

% Read the control file and generate stimulus trigger timestamps
switch lower(ext)
    case '.ncs'
        [ctrlTS, dg_Nlx2Mat_Samples] = dg_readCSC( ...
            ctrlfilepath );
    case '.mat'
        load(ctrlfilepath);
        ctrlTS = dg_Nlx2Mat_Timestamps; %#ok<NODEF>
    otherwise
        error('dg_rmStimArtifact:badinfiletype', ...
            'Unrecognized file extension: %s', ext);
end
if verboseflag
    fprintf('Read %s\n', ctrlfilepath);
end
sampleperiod = median(diff(ctrlTS)) / ...
    size(dg_Nlx2Mat_Samples,1);
allTS = round(repmat(ctrlTS, size(dg_Nlx2Mat_Samples,1), 1) ...
    + sampleperiod * repmat( ...
    (0:size(dg_Nlx2Mat_Samples,1) - 1)', 1, size(dg_Nlx2Mat_Samples, 2) ));
if isempty(threshvalopt)
    [threshval, isgood] = ctrlfileStats(dg_Nlx2Mat_Samples, thresh);
    if ~isgood
        error('dg_rmStimArtifact:ctrl', ...
            'The control file appears not to contain a control signal:\n%s', ...
            ctrlfilepath);
    end
else
    threshval = threshvalopt;
end
% trigsamp is first strictly-above-threshold sample
trigsamp = find( (dg_Nlx2Mat_Samples(2 : end) > threshval) ...
    & (dg_Nlx2Mat_Samples(1 : end - 1) <= threshval) ) + 1;
timeoutpts = round(timeout/sampleperiod);
trigsampdiff = diff(trigsamp);
trigsamp(find(trigsampdiff < timeoutpts) + 1) = [];
trigTS = reshape(allTS(trigsamp), [], 1); % TS of trigsamp
if verboseflag
    fprintf('Found %d trigger times\n', length(trigTS));
end

% Interpolate away the stimulus artifacts
sameTS = false(size(filenames));
for fileidx = 1:length(filenames)
    [~, ~, ext] = fileparts(filenames{fileidx});
    switch lower(ext)
        case '.ncs'
            [dg_Nlx2Mat_Timestamps, dg_Nlx2Mat_Samples, header] = ...
                dg_readCSC(fullfile(src, filenames{fileidx}));
            if vscale
                for k = 1:length(header)
                    if regexp(header{k}, '^\s*-ADBitVolts\s+')
                        ADBitVoltstr = regexprep(header{k}, ...
                            '^\s*-ADBitVolts\s+', '');
                        ADBitVolts = str2double(ADBitVoltstr);
                        if isempty(ADBitVolts)
                            warning('dg_rmStimArtifact:badADBitVolts', ...
                                'Could not convert number from:\n%s', ...
                                header{k} );
                        else
                            dg_Nlx2Mat_Samples = ADBitVolts ...
                                * dg_Nlx2Mat_Samples;
                            dg_Nlx2Mat_SamplesUnits = 'V'; %#ok<NASGU>
                            threshval = ADBitVolts * threshval;
                        end
                    end
                end
            else
                dg_Nlx2Mat_SamplesUnits = 'AD'; %#ok<NASGU>
            end
        case '.mat'
            load(fullfile(src, filenames{fileidx}));
        otherwise
            error('dg_rmStimArtifact:badinfiletype', ...
                'Unrecognized file extension: %s', ext);
    end
    if verboseflag
        fprintf('Read %s\n', fullfile(src, filenames{fileidx}));
    end
    if isempty(threshvalopt)
        [prethreshval, preisgood] = ctrlfileStats(dg_Nlx2Mat_Samples, thresh);
        if ~preisgood
            warning('dg_rmStimArtifact:noartifact', ...
                '%s may not actually contain stim artifacts', ...
                filenames{fileidx});
        end
    end
    sameTS(fileidx) = isequal(dg_Nlx2Mat_Timestamps, ctrlTS);
    sampleperiod = median(diff(dg_Nlx2Mat_Timestamps)) / ...
        size(dg_Nlx2Mat_Samples,1);
    numIntrp = round(diff(offsets) / sampleperiod) + 1;
    allTS = repmat(dg_Nlx2Mat_Timestamps, size(dg_Nlx2Mat_Samples,1), 1) ...
        + sampleperiod * repmat((0:size(dg_Nlx2Mat_Samples,1) - 1)', 1, ...
        size(dg_Nlx2Mat_Samples, 2));
    laststimsample = 1;
    for stimnum = 1:length(trigTS)
        startsample = laststimsample;
        endsample = min( numel(dg_Nlx2Mat_Samples), ...
            startsample + floor(5e6/sampleperiod) );
        if sameTS(fileidx)
            trigidx = trigsamp(stimnum);
        else
            trigidx = dg_binsearch(allTS(:), trigTS(stimnum), ...
                startsample, endsample);
        end
        if allTS(trigidx) - trigTS(stimnum) > sampleperiod/2
            index = trigidx - 1;
        else
            index = trigidx;
        end
        index = index + round(offsets(1) / sampleperiod);
        if index < 1
            warning('dg_rmStimArtifact:tooearly', ...
                'The artifact at %d s is too early to interpolate.', ...
                trigTS(stimnum)*1e-6 );
            continue
        end
        if index > numel(dg_Nlx2Mat_Samples)
            warning('dg_rmStimArtifact:toolate', ...
                'The artifact at %d s is too late to interpolate.', ...
                trigTS(stimnum)*1e-6 );
            continue
        end
        if index < numel(allTS)+1
            if index+numIntrp < numel(dg_Nlx2Mat_Samples)+1
                dg_Nlx2Mat_Samples(index + (0:numIntrp)) = linspace( ...
                    dg_Nlx2Mat_Samples(index), ...
                    dg_Nlx2Mat_Samples(index+numIntrp), numIntrp+1 );
            else
                dg_Nlx2Mat_Samples(index:end) = dg_Nlx2Mat_Samples(index);
            end
        else
            break
        end
        laststimsample = index;
    end
    if isempty(threshvalopt)
        maxfraclarge = max([1e-6 10/numel(dg_Nlx2Mat_Samples)]);
        if preisgood && sum(dg_Nlx2Mat_Samples(:) > prethreshval) >= ...
                numel(dg_Nlx2Mat_Samples) * maxfraclarge
            warning('dg_rmStimArtifact:artifact', ...
                'At least %5d of the samples in %s have large values after interpolation.', ...
                maxfraclarge, filenames{fileidx});
        end
    end
    [~, n] = fileparts(filenames{fileidx});
    matfilepath = fullfile(dest, [n, 'rmstim.mat']);
    save(matfilepath, 'dg_Nlx2Mat_Timestamps', ...
        'dg_Nlx2Mat_Samples', 'dg_Nlx2Mat_SamplesUnits', '-v7.3');
    if verboseflag
        fprintf('Wrote %s\n', matfilepath);
    end
end

end


function [threshval, isgood] = ctrlfileStats(samples, thresh)
%OVERVIEW
% Recognizes and returns values from files that meet one of the following
% descriptions:
%   1. Sample value histogram has its largest peak near (within one tenth
% of the maximum excursion from) the median, and at least <fracclipped> of
% the samples are at the maximum sample value, indicating that they are
% most likely hard-clipped and the waveform can be treated as a rectangular
% waveform;
%   2. Sample value histogram has its largest peak as stated above, and a
% another peak at higher sample values than the first, with at least
% <peakthresh> as many counts as the largest peak, indicating that the
% there is a "high" state that is soft-clipped (e.g. as a result of lowpass
% filtering), and at least <fracclipped> of the samples are above <medval>
% + <thresh> * (<highval> - <medval>), where <highval> is the sample value
% at the second peak;
%   3. Sample value histogram has largest peak as stated above, and at
% least <fracclipped> of the samples are more than <ratio> * <range from
% 10th prctile to 90th prctile> above <medval>, suggesting that although it
% is not clipped, there is some sort of large transient signal that could
% be used as a control signal.
%INPUTS
% samples: the data.
% thresh: same as <thresh> in dg_rmStimArtifact; see OVERVIEW.
% peakthresh: the fraction of the bin count at the first mode that is
%   acceptable for the "second mode" (if there is one).
%OUTPUTS
% threshval:
%   Case 1. <medval> + <thresh> * (<maxval> - <medval>).
%   Case 2: <medval> + <thresh> * (<highval> - <medval>), where <highval>
%       is the sample value at the second largest peak.
%   Case 3. <medval> + <thresh> * (<maxval> - <medval>).
%   Otherwise: NaN.
% isgood: true if the stats fit one of the cases described in OVERVIEW;
%   false otherwise.
isgood = false;
threshval = NaN;
nbins = 100;
peakthresh = 1e-4;
fracclipped = max([1e-6 10/numel(samples)]);
medval = median(samples(:));
maxval = max(samples(:));
minval = min(samples(:));
valrange = maxval - minval;
binedges = minval + (0:nbins) * valrange / nbins;
binctrs = binedges + 0.5 * valrange / nbins;
counts = histc(samples(:), binedges);
[peakcnt, maxidx] = max(counts);
modeval = binctrs(maxidx);
if modeval - medval > medval + (maxval - medval) / 10
    % The mode is not near the median
    return
end
% <modeval> satisfies the largest peak criteria.


if sum(samples(:) == maxval) / numel(samples) >= fracclipped
    % Case 1.
    threshval = medval + thresh * (maxval - medval);
else
    % Check for Case 3.  In a Gaussian distribution, the 1e-6 point is
    % reached at 4.75 SD from the mean, whereas the 10% point is at 1.28 SD
    % from the mean; so if the size of the tail above (4.75 / 1.28) times
    % the 90th prctile minus the median is greater than 1e-6, then the
    % signal is non-Gaussian in the way that a control signal should be.
    % We throw in a factor 2 as a safety margin.
    ratio = 4.75 / 1.28;
    if sum( samples(:) > 2 * ratio * (prctile(samples(:), 90) - medval) ...
            + medval ) / numel(samples) >= 1e-6
        % Case 3 confirmed.
        threshval = medval + thresh * (maxval - medval);
    else
        % Check for Case 2.  This requires that we smooth the histo.
        sw = hanning(7);
        counts = conv(counts, sw, 'same') / sum(sw);
        % Find the <startidx>, the bin index at the valley
        % between modes.  We assume that the non-control component of the
        % signal (LFP, noise, whatever) has a Gaussian distribution as for
        % Case 3:
        exclusionval = ratio * (prctile(samples(:), 90) - medval);
        [~, startidx] = min(abs(binctrs - exclusionval));
        while startidx < length(counts) && counts(startidx) ...
                > counts(startidx+1)
            startidx = startidx + 1;
        end
        if startidx < length(counts)
            % Found valley.
            [peakcnt2, maxidx2] = max(counts(startidx:end));
            maxidx2 = maxidx2 + startidx - 1;
            if peakcnt2 >= peakthresh * peakcnt
                % second peak is big enough
                case2threshval = medval + ...
                    thresh * (binctrs(maxidx2) - medval);
                if sum(samples(:) > case2threshval) / numel(samples) ...
                        >= fracclipped
                    % Case 2 confirmed.
                    threshval = case2threshval;
                else
                    return
                end
            else
                return
            end
        else
            % There is no valley.
            return
        end
    end
end
isgood = true;
    
end



