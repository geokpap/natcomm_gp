function [freqs, sumP, sumsqrsP, N, matrix, Wraw] = dg_crossFreqMorlet(...
    windows, rawsamples, varargin)
% Applies Morlet-based wavelet transform to each column of <windows>, then
% computes correlations between pseudo-power at different frequencies
% across the columns of <windows>, i.e. using each sample time as a
% separate observation.  If the return value <matrix> is not used, it is
% not computed.
%INPUTS
% windows: 2-column array of sample indices as returned by
%   lfp_collectDataWindows.
% rawsamples: cell vector of sample windows.  Each cell must contain a
%   numeric array of the same size.
%OUTPUTS
% freqs: column vector of pseudo-frequency for each row and column of
%     <matrix>. If 'phase' option is used, a second column contains the
%     phase values for each column of <matrix>.
% sumP: column vector of the sum across trials of the magnitude of the
%     wavelet transform. Contains two columns if two <filenums> are
%     given (see lfp_BLspectrum for additional notes).
% sumsqrsP: column vector of the sum across trials of the squared
%     magnitude of the wavelet transform. Contains two columns if two
%     <filenums> are given.
% N: the number of trials that went into <sumP> and <sumsqrsP>.
% matrix: cross-covariance matrix between channels of the magnitudes
%     for each value of <f> (auto-covariance if only one channel is given).
%     Rows refer to <f> in rawsamples{1}, columns refer to <f> in
%     rawsamples{2}.
% Wraw: cell row vector of the complex morletgrams for each window in
%   the last batch, in {window}(freq, timepoint) format.  
%OPTIONS
% 'detrend', detrendflag
%   If <detrendflag> is true (default), applies Matlab 'detrend' function
%   to remove linear trend separately from each column of <windows> before
%   processing.  Note that this implicitly includes the equivalent of
%   'rmdc' as part of the detrending.  If false, there is no detrending or
%   removal of dc.  If <detrendflag> is 'rmdc', removes dc without removing
%   trend.
% 'freqlim', freqlim
%   <freqlim> is a 2-element row vector specifying the lower and upper
%   limits of frequencies to include in the correlation matrix.  These are
%   specified in the same units as f (see 'fs' option below).  If not
%   given, [1 fs/2] is used.
% 'fs', samplefreq
%   Converts frequencies to Hertz based on <samplefreq>.  If not given, fs
%   is assumed to be 1.
% 'fspacing', fspacing
%   ratio between scales on successive rows, which gets passed in to
%   dg_morletgram.  Use a value of 10^(1/n) where <n> is an integer if you
%   want nice round frequency tick mark labels. Default: 10^(1/32).
% 'param', value
%   the Morlet wavelet width parameter. This is approximately the number of
%   cycles in the wavelet. Default = 6.
% 'phase'
%   use the phase instead of the power of the wavelet transforms from
%   channel 2.
% 'pts2trim', numpts - in cases where the <windows> have had some margin
%   time added before and after, it is necessary to trim the margin off of
%   the wavelet transforms before running the statistics.  <numpts> is the
%   number of points to remove at each end.
% 'verbose'
%   optional.
%NOTES
% When computation must be done in multiple batches, intermediate results
% are saved in the current working directory in temporary files with names
% starting with 'dg_crossFreqCorr_temp'.  If an error is encountered before
% the temporary files are deleted, they must be deleted manually.
%   STATISTICS: since every sample is treated as a separate observation,
% there is a built-in autocorrelation that has a pink character, i.e. it
% gets monotonically smaller with pseudo-frequency.  This means that
% Gaussian white noise should have a pink auto-crossFreqMorlet correlation
% spectrum.

%$Rev: 257 $
%$Date: 2017-06-12 14:41:30 -0400 (Mon, 12 Jun 2017) $
%$Author: dgibson $

detrendflag = true;
freqlim = [];
fs = 1;
fspacing = 10^(1/32);
verboseflag = false;
param = 6;
phaseflag = false; %#ok<NASGU>
pts2trim = 0;
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'detrend'
            argnum = argnum + 1;
            detrendflag = varargin{argnum};
        case 'freqlim'
            argnum = argnum + 1;
            freqlim = varargin{argnum};
        case 'fs'
            argnum = argnum + 1;
            fs = varargin{argnum};
        case 'param'
            argnum = argnum + 1;
            param = varargin{argnum};
        case 'fspacing'
            argnum = argnum + 1;
            fspacing = varargin{argnum};
            if fspacing <= 1
                error('dg_crossFreqMorlet:fspacing', ...
                '''fspacing'' must be greater than 1.');
            end
        case 'phase'
            phaseflag = true; %#ok<NASGU>
            error('dg_crossFreqMorlet:phase', ...
                'The ''phase'' option is not implemented.');
        case 'pts2trim'
            argnum = argnum + 1;
            pts2trim = varargin{argnum};
        case 'verbose'
            verboseflag = true;
        otherwise
            error('dg_crossFreqMorlet:badoption', ...
                'The option %s is not recognized.', ...
                dg_thing2str(varargin{argnum}));
    end
    argnum = argnum + 1;
end

if isempty(freqlim)
    freqlim = [1 fs/2];
end
if freqlim(1) == 0
    error('dg_crossFreqMorlet:freq', ...
        'The low frequency limit must be nonzero.');
end
if freqlim(1) >= freqlim(2)
    error('dg_crossFreqMorlet:freq2', ...
        'The low frequency limit must be lower than the high frequency limit.');
end

numch = length(rawsamples);
if numch == 0
    error('dg_crossFreqMorlet:numch', ...
        '<rawsamples> is empty.');
elseif numch > 2
    error('dg_crossFreqMorlet:rawsamples', ...
        '<rawsamples> contains too  many channels.');
end
if any(windows(:) > numel(rawsamples{1}))
    error('dg_crossFreqMorlet:windows', ...
        '<windows> exceeds bounds of <rawsamples>.');
end
for ch = 2:numch
    if ~isequal(size(rawsamples{1}), size(rawsamples{ch}))
        error('dg_crossFreqMorlet:rawsamples2', ...
            '<rawsamples> contains channels of different sizes.');
    end
    if any(windows(:) > numel(rawsamples{ch}))
        error('dg_crossFreqMorlet:windows2', ...
            '<windows> exceeds bounds of <rawsamples>.');
    end
end

% convert everything to log units to construct the frequency grid.
logfspacing = log(fspacing);
logfreqlim = log(freqlim);
logf = logfreqlim(1) : logfspacing : logfreqlim(2);
freqs = exp(logf);
numfreqs = length(freqs);
FourierFactor = 4*pi/(param+sqrt(2+param^2));
scales = 1./(FourierFactor * freqs);
dt = 1/fs;

% The computation is a two-pass process, and is batched that way.  The
% first pass consists of computing wavelet transforms and pseudopower
% variances using the lfp_BLspectrum method of accumulating sumP and
% sum-of-squares.  If more than one batch is needed, each batch of wavelet
% transforms is saved to a temp file.  The second pass consists of
% subtracting the mean wavelet transform from each batch of spectra and
% accumulating the cross-covariance sum batch by batch.
%

% Batch according to largest available contiguous block of memory (adapted
% from dg_downsample)
tempfilenum = 1;
tempbase = sprintf('dg_crossFreqMorlet_temp%d_', tempfilenum);
while ~isempty(dir([tempbase '*']))
    tempfilenum = tempfilenum + 1;
    tempbase = sprintf('dg_crossFreqMorlet_temp%d_', tempfilenum);
end
if isunix
    % The 'memory' function does not exist on Mac, so just assume we can
    % use one half of the maximum possible after subtracting easily
    % determinable usage, which by definition is the usage calculated by
    % this snippet adapted from the Matlab support website, operating in
    % the current workspace and in the base workspace.  Unfortunately, I
    % have no idea how to find out how much memory is actually available on
    % Mac/FreeBSD, but
    % http://developer.apple.com/documentation/Performance/Conceptual/Manag
    % ingMemory/Articles/AboutMemory.html says that each process can have
    % up to 4 GB of addressable space; of course 1.5 GB of that is consumed
    % by Matlab overhead, leaving, oh, let's just call the maximum possible
    % "2 GB", i.e. 2^31. 
    mem_elements = whos;
    if size(mem_elements,1) > 0
        for i = 1:size(mem_elements,1)
            memory_array(i) = mem_elements(i).bytes; %#ok<*AGROW>
        end
        memory_in_use = sum(memory_array);
    else
        memory_in_use = 0;
    end
    mem_elements = evalin('base','whos');
    if size(mem_elements,1) > 0
        for i = 1:size(mem_elements,1)
            memory_array(i) = mem_elements(i).bytes;
        end
        
        memory_in_use = memory_in_use + sum(memory_array);
    end
    bytesinblock = (2^31 - memory_in_use) / 2;
else
    % It must be Windows
    userview = memory;
    bytesinblock = userview.MaxPossibleArrayBytes;
end
% I've included a fudge factor of 2 in case Matlab feels
% a need to make a copy of everything:
maxbytes2use = bytesinblock / 2;
if bytesinblock < 24
    error('dg_crossFreqMorlet:outofmem', ...
        'There are only %d bytes available; pack memory', bytesinblock);
end
% calculate the maximum size of each batch (i.e. the number of samples in
% the batch) and the estimated number of batches <estbatches>.  I assume
% that exactly 2 waveforms worth of memory are needed to store individual
% waveforms during the computations.
winbytespersample = (length(freqs) + 2) * 8;
batchsize = fix(maxbytes2use/winbytespersample);
if numch == 2
    % Two-channel version
    batchsize = fix(batchsize/2);
end

%
% First pass
%

totsamples = sum(diff(windows, 1, 2));
ismultibatch = totsamples > batchsize;
% We cannot know a priori the number of windows in each batch because
% windows are of variable size, nor the exact number of batches.  But we
% can estimate the number of batches if we want to be verbose.
if verboseflag
    estnumbatches = ceil(totsamples/batchsize);
    fprintf('This will take approx. %d batches\n', estnumbatches);
end
% address as sumP(freqidx, chidx); sumsqrsP likewise.
sumP = zeros(numfreqs, numch);
sumsqrsP = zeros(numfreqs, numch);
N = 0; % total number of samples processed per channel
numwins = size(windows,1);
tempfilenames = {};
batchnum = 1;
W = {}; % accumulates wavelet transforms within a batch
startwin = 1;
endwin = startwin;
batchwinidx = 0;
numsampinbatch = 0;
if verboseflag
    startT = tic;
end
for winidx = 1:numwins
    % compute one window's wavelet transform
    % in W{batchwinidx}(freqs, samples, ch).
    sampidx = windows(winidx,1):windows(winidx,2); % into <rawsamples>
    winlen = length(sampidx);
    numgoodpts = winlen - 2 * pts2trim;
    isgoodpt = [false(1, pts2trim) true(1, numgoodpts) false(1, pts2trim)];
    numsampinbatch = numsampinbatch + numgoodpts;
    N = N + numgoodpts;
    batchwinidx = batchwinidx + 1;
    % Force W to be a row vector.
    W{1, batchwinidx} = zeros(numfreqs, numgoodpts, numch);
    for chidx = 1:numch
        if isequal(detrendflag, true)
            sig1 =  struct( ...
                'val', detrend(rawsamples{chidx}(sampidx)), ...
                'period', dt );
        elseif isequal(detrendflag, 'rmdc')
            sig1 =  struct( ...
                'val', rawsamples{chidx}(sampidx) - ...
                nanmean(rawsamples{chidx}(sampidx)), ...
                'period', dt );
        else
            sig1 =  struct( ...
                'val', rawsamples{chidx}(sampidx), ...
                'period', dt );
        end
        wt = cwtft(sig1, 'scales', scales, 'wavelet', {'morl', param});
        W{1, batchwinidx}(:, :, chidx) = wt.cfs(:, isgoodpt);
        % sum the pseudopower over all timepoints in each window:
        sumP(:, chidx) = sumP(:, chidx) + ...
            sum(wt.cfs(:, isgoodpt) .* conj(wt.cfs(:, isgoodpt)), 2);
        sumsqrsP(:, chidx) = sumsqrsP(:, chidx) + ...
            sum((wt.cfs(:, isgoodpt) .* conj(wt.cfs(:, isgoodpt))) .^ 2, 2);
    end
    batchdone = (winidx == numwins) || ...
        numsampinbatch + diff(windows(winidx+1, :)) + 1 > batchsize;
    if batchdone
        if nargout >= 6
            Wraw = W;
        end
        W = cell2mat(W); % a desperate attempt to conserve memory
        if ismultibatch
            % Save this batch and initialize for next one.
            tempfilenames{end+1} = sprintf( '%s%07d.mat', tempbase, ...
                startwin );
            % concatenate the cell contents using cell2mat, leaving just
            % W(freqs, samples, ch).
            if nargout >= 5
                % Save temp files for computing <matrix>.
                save(tempfilenames{end}, 'W');
                if verboseflag
                    fprintf('Saved file %s\n', tempfilenames{end} );
                end
            end
            % initialize for next batch:
            batchnum = batchnum + 1;
            W = {};
            startwin = endwin + 1;
            endwin = startwin;
            batchwinidx = 0;
            numsampinbatch = 0;
        end
        if verboseflag
            fprintf( ...
                'Batch %.0f done, ending window %.0f\n', ...
                batchnum - 1, startwin - 1 );
            T = toc(startT);
            Trem = (T / (batchnum - 1)) * (estnumbatches - batchnum + 1);
            fprintf('Est. time remaining: %.1f minutes\n', Trem/60);
        end
    end
    endwin = endwin + 1;
end

% avgP(freqidx, chidx) is like sumP:
avgP = sumP/N;
% We use the sqrt-sum-of-squared-diffs, i.e. sqrt(N-1)*std, instead of
% std for sigma, which must be what Masimore ... Redish 2004 eqn 1
% meant, because otherwise the final value is proportional to N-1. 
%   (This in fact turns out to be correct; see corrigendum in J Neurosci
%   Methods. 2005 Feb 15;141(2):333. - DG 16-Sep-2015)
% Using varP = (sumsqrsP - (sumP).^2/N)/(N-1), we get the formula for
% sigma; address as sig(freqidx, chidx):
sig = sqrt(sumsqrsP - sumP.^2/N);
% idiot checquerei on the above, starting with the formula from
% https://en.wikipedia.org/wiki/Standard_deviation#Rapid_calculation_methods:
% s = sqrt((N*s2 - s1^2)/(N*(N-1)));
% s = sqrt((N*sumsqrsP - sumP^2)/(N*(N-1)));
% s = sqrt((sumsqrsP - sumP^2/N)/(N-1));
% s = sqrt(1/(N-1)) * sqrt(sumsqrsP - sumP^2/N);
% sqrt(N-1)*s = sqrt(sumsqrsP - sumP^2/N);
% Checque (please).

if nargout < 5
    % <matrix> is not used, so don't compute it.
    return
end

%
% Second pass:
% accumulating the cross-covariance sum <xcovsum> batch by
% batch.  For each batch, the contribution to xcovsum(I,J) is the dot
% product of the time series of power re: avg at freq I in the first
% channel with the time series of power at freg J in the second.
%

xcovsum = zeros(numfreqs);
if ismultibatch
    numbatches = length(tempfilenames);
else
    numbatches = 1;
end
for tempfilenum = 1:numbatches
    if ismultibatch
        if verboseflag
            fprintf( ...
                'Loading temp file #%.0f: %s\n', ...
                tempfilenum, tempfilenames{tempfilenum} );
        end
        load(tempfilenames{tempfilenum});
        delete(tempfilenames{tempfilenum});
    end
    % address as W(freqs, samples, ch).
    % address as avgP(freqs, ch), diffs(freqs, samples, ch).
    diffs = W .* conj(W) - repmat( ...
        reshape(avgP, [size(avgP,1) 1 size(avgP,2)]), ...
        1, size(W,2) );
    if numch == 2
        % Two channel xcov
        xcovsum = xcovsum + diffs(:,:,1) * diffs(:,:,2)';
    else
        % One channel, i.e. auto-xcov
        xcovsum = xcovsum + diffs * diffs';
    end
end
% <xcovsum> now contains the numerator of eqn 1, summed over all batches.
sigmatrix = NaN(length(sig));
for I = 1:length(sig)
    for J = 1:length(sig)
        if numch == 2
            sigmatrix(I,J) = sig(I, 1) * sig(J, 2);
        else
            sigmatrix(I,J) = sig(I) * sig(J);
        end
    end
end
matrix = xcovsum ./ sigmatrix;

