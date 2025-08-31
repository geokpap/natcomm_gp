function [corr, f, CV, hI, hCB] = dg_crossFreqCorr(data, varargin)
%[corr, f] = dg_crossFreqCorr(data)
% Applies multi-taper spectral analysis to each column of <data>, then
% computes correlations between spectral power at different frequencies
% across the columns of <data>.  Creates a figure and plots the
% cross-frequency cross-covariance matrix in pseudocolor, with frequency
% index for channel 1 on the Y axis and freq idx of channel 2 on the X.
%INPUTS
% data: contains one time window of consecutive samples in each column.
%   The numbers of rows (samples per window) and columns (windows) are both
%   variable.
%OUTPUTS
% corr: actually the crosscovariance, which ranges from -1 to 1.
% f: frequency for each row and column of <corr>.
% CV: coefficient of variation for each frequency in <f>.
%   When 'getfrom'/'getfrom2' are being used, CV has two columns, one for
%   each channel.
% hI, hCB: from dg_showGram.
%OPTIONS
%'descrip', descrip
%  <descrip> is a string that is included as an additional line in the
%  figure title.
%'detrend', detrendflag
%  If <detrendflag> is true (default), applies Matlab 'detrend' function to
%  remove linear trend separately from each column of <data> before
%  processing.  Note that this implicitly includes the equivalent of 'rmdc'
%  as part of the detrending.  If false, there is no detrending or removal
%  of dc.  If <detrendflag> is 'rmdc', removes dc without removing trend.
%'freqlim', freqlim
%  <freqlim> is a 2-element row vector specifying the lower and upper
%  limits of frequencies to include in the correlation matrix.  These are
%  specified in the same units as f (see 'fs' option below).
%'fs', samplefreq
%  Converts frequencies to Hertz based on <samplefreq>.  If not given,
%  f = 0:(nfft-1) is returned, where <nfft> is the number of data points
%  per window after padding.
%'getfrom', datasource
%  <datasource> is a large and unwieldy time series, portions of which are
%  analyzed in batches, and we do NOT want Matlab to copy it.  In this
%  case, <data> does NOT contain the actual data, but is a 2-column array
%  that contains starting point indices in the first column and ending
%  point indices in the second column, constrained s.t. the number of
%  points between the starting point and ending point is always the same.
%  The indices point into <datasource>. Neither <datasource> nor <data> are
%  altered in any way, so Matlab "should" pass them by reference (i.e., not
%  copy them).
%'getfrom2', datasource
%  Same as 'getfrom', but supplies data for a second channel and invokes
%  cross-channel cross-frequency correlation computation.
%'k', numtapers
%  Forces number of multitapers to <numtapers>; default is 1.
%'noplot'
%  Supresses plotting of correlation matrix.
%'nw', nw
%  Forces "time-bandwidth product" setting for multitapers to <nw>; default
%  is 1.8.
%'pad', N
%  N=0 does not pad at all (note that this differs from lfp_mtspectrum),
%  N=1 pads to twice the next greater power of 2, N=2 pads to four times
%  that length, etc.  The default value is N=0.
% 'verbose' - display too much running status information
%NOTES
% When computation must be done in multiple batches, intermediate results
% are saved in the current working directory in temporary files with names
% starting with 'dg_crossFreqCorr_temp'.  If an error is encountered before
% the temporary files are deleted, they must be deleted manually.
%
% This is copied from Rev 36, which is the one that actually works as of
% 16-Dec-2010.
%REFERENCES
% Masimore ... Redish 2004: J Neurosci Methods. 2004 Sep 30;138(1-2):97-105.

%$Rev: 226 $
%$Date: 2015-09-17 19:45:40 -0400 (Thu, 17 Sep 2015) $
%$Author: dgibson $

datasource = [];
datasource2 = [];
descrip = '';
detrendflag = true;
freqlim = [];
fs = 0;
K = 1;
NW = 1.8;
numch = 1;
pad = 0;
plotflag = true;
verboseflag = false;
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'descrip'
            argnum = argnum + 1;
            descrip = varargin{argnum};
        case 'detrend'
            argnum = argnum + 1;
            detrendflag = varargin{argnum};
        case 'freqlim'
            argnum = argnum + 1;
            freqlim = varargin{argnum};
        case 'fs'
            argnum = argnum + 1;
            fs = varargin{argnum};
        case {'getfrom' 'getfrom2'}
            optstr = varargin{argnum};
            argnum = argnum + 1;
            if isempty(varargin{argnum})
                error('dg_crossFreqCorr:datasource', ...
                    '<datasource> is empty');
            end
            if isequal(optstr, 'getfrom')
                datasource = varargin{argnum};
            else
                datasource2 = varargin{argnum};
                numch = 2;
            end
            if isempty(data)
                error('dg_crossFreqCorr:data', ...
                    '<data> is empty');
            end
            if ~all(data(:,2) - data(:,1) == data(1,2) - data(1,1))
                error('dg_crossFreqCorr:data2', ...
                    '<data> windows are not all the same length');
            end
        case 'k'
            argnum = argnum + 1;
            K = varargin{argnum};
        case 'noplot'
            plotflag = false;
        case 'nw'
            argnum = argnum + 1;
            NW = varargin{argnum};
        case 'pad'
            argnum = argnum + 1;
            if ~isa(varargin{argnum}, 'double')
                error('lfp_spec:badpad', ...
                    'The padding factor must be a number');
            end
            pad = varargin{argnum};
        case 'verbose'
            verboseflag = true;
        otherwise
            error('dg_crossFreqCorr:badoption', ...
                ['The option "' varargin{argnum} ...
                '" is not recognized.'] );
    end
    argnum = argnum + 1;
end

% "Time-bandwidth product NW must be less than N/2":
if NW >= size(data,1) / 2
    error('dg_crossFreqCorr:toofewpts', ...
        'There are too few rows in <data> for the specified value of nw');
end

if isequal(detrendflag, true)
    detrendstr = 'y';
elseif isequal(detrendflag, 'rmdc')
    detrendstr = 'rmdc';
else
    detrendstr = 'no';
end

if isempty(datasource)
    % not using 'getfrom'/'getfrom2'
    numwins = size(data, 2);
    winsamples = size(data, 1);
else
    numwins = size(data, 1);
    winsamples = data(1,2) - data(1,1) + 1;
end
tapers = dpss(winsamples, NW, K);
if pad > 0
    nfft = 2^(nextpow2(winsamples) + pad);
else
    nfft = winsamples;
end
f = 0:(nfft-1);
if fs ~= 0
    f = fs * f / nfft;
end
if ~isempty(freqlim)
    if ~isequal(size(freqlim), [1 2])
        error('dg_crossFreqCorr:freqlim', ...
            '<freqlim> must be a two-element row vector.');
    end
    isfreq2use = (f >= freqlim(1) & f <= freqlim(2));
    f = f(isfreq2use);
else
    isfreq2use = true(size(f));
end
numfreqs2use = sum(isfreq2use);

%
% Using Masimore ... Redish 2004 eqn 1, where <corr> is their rho and <P>
% is their S, k corresponds to the columns of P (but using
% sqrt-sum-of-squared-diffs instead of std for sigma, which must be what
% they meant, because otherwise the final value is proportional to N-1).
%
% The computation is a two-pass process, and is batched that way.  The
% first pass consists of computing power spectra and variances using the
% lfp_BLspectrum method of accumulating sum and sum-of-squares.  If more
% than one batch is needed, each batch of power spectra is saved to a temp
% file.  The second pass consists of subtracting the mean spectrum from
% each batch of spectra and accumulating the cross-covariance sum batch by
% batch.
%

% Batch according to largest available contiguous block of memory (adapted
% from dg_downsample)
tempfilenum = 1;
tempbase = sprintf('dg_crossFreqCorr_temp%d_', tempfilenum);
while ~isempty(dir([tempbase '*']))
    tempfilenum = tempfilenum + 1;
    tempbase = sprintf('dg_crossFreqCorr_temp%d_', tempfilenum);
end
winbytes = winsamples * 8;
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
startwin = 1;   % first window in current batch
if bytesinblock < 24
    error('dg_crossFreqCorr:outofmem', ...
        'There are only %d bytes available; pack memory', bytesinblock);
end
% calculate estimated number of batches <estbatches>.

% Some book on Google (ISBN 9971509199) claims on p. 333 that "a direct
% DFT requires N^2 complex multiplications", so it can't possibly
% require more storage than O(N^2). Let's assume it's O(N), which would
% make life simple.
%   The data themselves require <winbytes>, and the indices into
% <datasource> require the same, so that's 2*winbytes required per
% window.  I seem to recall Matlab feels compelled to copy data when
% indices are applied, so that makes 3*winbytes.  The method of
% calculating the indices involves adding two arrays of the same size,
% bringing us up to 5*winbytes. Then there are the tapered windows
% (data_proj), so that's (K+5) * winbytes. Aside from producing a final
% power spectrum, we also need working storage to do an FFT on <K>
% tapers per window and average them, so that makes (2*K+6) * winbytes.
% We also need to store the tapers themselves, so that makes (3*K+6) *
% winbytes.  And it turns out that the "permute" in chronux_mtfftc can cause
% problems, so that needs another copy of data, but the
% "data=data(:,:,ones(1,K));" multiplies the size of data by K, so what
% we finally need in the end for the "permute" is 2*K*winbytes,
% bringing the grand total to (5*K+6) * winbytes.
%   The "P = squeeze..." requires twice the size of S to assign
% (liberally considering P to be the same size as S, which it isn't
% really), three times counting the "conj", four times counting
% temporary storage for the product, and another K*nfft*8 for the "mean".
% S is K * nfft * 8 bytes per window, totalling 4*batchsize*K*nfft*8 +
% K*nfft bytes required.  This means we must compute batchsize twice,
% and use the smaller value.
batchsize = fix(min( ...
    maxbytes2use/((5*K+6) * winbytes), ...
    (maxbytes2use - K*nfft*8)/(4*K*nfft*8) ));
if ~isempty(datasource2)
    % Two-channel version
    batchsize = fix(batchsize/2);
end
endwin = batchsize; % last window in current batch
% Unlike dg_downsample, there is no overlap of batches.  However, out
% of laziness, I retain the same caveat, which is probably no longer true:
% "Since <estbatches> is only approximate, it can NOT be used as a loop
% termination criterion."
estbatches = 2 + ceil((numwins - 2*(batchsize-1))/(batchsize-2));
ismultibatch = estbatches > 1;
if ~ismultibatch
    endwin = numwins;
end
if verboseflag
    fprintf('Processing in approx. %d batches\n', estbatches);
end

% First pass
done = false;
% address as sumP(freqidx, chidx); sumsqrsP likewise.
sumP = zeros(numfreqs2use, numch);
sumsqrsP = zeros(numfreqs2use, numch);
N = 0;
batchnum = 0;
tempfilenames = {};
while ~done
    batchnum = batchnum + 1;
    if verboseflag
        fprintf( ...
            'Batch %.0f, windows %.0f-%.0f\n', ...
            batchnum, startwin, endwin );
    end
    % compute one batch of power spectra
    numwin = endwin - startwin + 1;
    if isempty(datasource)
        batchsamples = data(:, startwin:endwin);
    else
        % Using 'getfrom'/'getfrom2':
        samplidx = repmat((0:winsamples-1)', 1, numwin);
        startidx = repmat(data(startwin:endwin,1)', winsamples, 1);
        batchsamples = datasource(startidx + samplidx);
        if ~isempty(datasource2)
            batchsamples(:,:,2) = datasource2(startidx + samplidx);
        end
    end
    % <batchsamples> is now in windownum X samples X channel form.
    clear samplidx startidx
    if isequal(detrendflag, true)
        % Annoyingly, 'detrend' does not work on >2D arrays
        for k = 1:size(batchsamples,3)
            batchsamples(:,:,k) = detrend(batchsamples(:,:,k));
        end
    elseif isequal(detrendflag, 'rmdc')
        batchsamples = batchsamples - ...
            repmat(mean(batchsamples,1), size(batchsamples,1), 1);
    end
    % Annoyingly, 'chronux_mtfftc' does not work on >2D arrays.
    % S comes back in form frequency index x taper index x windows:
    S = chronux_mtfftc(batchsamples(:,:,1), tapers, nfft);
    if numch == 2
        S(:,:,:,2) = chronux_mtfftc(batchsamples(:,:,2), tapers, nfft);
    end
    clear batchsamples
    % Compute <P>  (spectral power), in [freq X window] format - i.e. each
    % column of <P> represents the same column of <data> (and each layer of
    % <P> represents the same channel of data)
    P = squeeze(mean(conj(S(isfreq2use,:,:,1)) .* S(isfreq2use,:,:,1), 2));
    if numch == 2
        P(:,:,2) = squeeze(mean(conj(S(isfreq2use,:,:,2)) .* S(isfreq2use,:,:,2), 2));
    end
    clear S
    % address as P(freqidx, windownum, channel), sumP(freqidx, chidx),
    % sumsqrsP(freqidx, chidx):
    sumP = sumP + squeeze(sum(P,2));
    sumsqrsP = sumsqrsP + squeeze(sum(P.^2,2));
    N = N + size(P, 2);
    if ismultibatch
        tempfilenames{end+1} = sprintf('%s%07d.mat', tempbase, startwin);
        if verboseflag
            fprintf( ...
                'Saving file %s\n', ...
                tempfilenames{end} );
        end
        save(tempfilenames{end}, 'P');
        if endwin == numwins
            done = true;
        else
            % update indices for next batch:
            startwin = endwin + 1;
            endwin = min(startwin + batchsize - 1, numwins);
        end
    else
        done = true;
    end
end
if N ~= numwins
    error('oops');
end
% like sumP, avgP(freqidx, chidx):
avgP = sumP/numwins;
% We use the sqrt-sum-of-squared-diffs, i.e. sqrt(N-1)*std, instead of
% std for sigma, which must be what Masimore ... Redish 2004 eqn 1
% meant, because otherwise the final value is proportional to N-1. 
%   (This in fact turns out to be correct; see corrigendum in J Neurosci
%   Methods. 2005 Feb 15;141(2):333. - DG 16-Sep-2015)
% Using varP = (sumsqrsP - (sumP).^2/N)/(N-1), we get the formula for
% sigma; address as sig(freqidx, chidx):
sig = sqrt(sumsqrsP - sumP.^2/numwins);
SD = sig/sqrt(N-1);
CV = SD ./ avgP;
    
% Second pass: accumulating the cross-covariance sum <xcovsum> batch by
% batch.  For each batch, the contribution to xcovsum(I,J) is the dot
% product of the time series of power re: avg at freq I in the first
% channel with the time series of power at freg J in the second.
xcovsum = zeros(numfreqs2use);
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
    % address as P(freqidx, windownum, channel), avgP(freqidx, chidx),
    % diffs(freqidx, windownum, channel):
    diffs = P - repmat(reshape(avgP, [size(avgP,1) 1 size(avgP,2)]), 1, size(P,2));
    % This is where the two-channel version becomes a headache. <diffs> is
    % in freq X window X channel form, and what we want to correlate is one
    % channel's diffs against the other one's:
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
corr = xcovsum ./ sigmatrix;

if plotflag
    hF = figure;
    if fs == 0
        winsizestr = sprintf('%.0f pts', winsamples);
    else
        winsizestr = sprintf('%g sec', winsamples / fs);
    end
    titlestr = sprintf( ...
        'Cross-Frequency Covariance, win=%s detrend=%s nw=%.2f k=%.0f, pad=%.0f', ...
        winsizestr, detrendstr, NW, K, pad);
    if ~isempty(descrip)
        titlestr = sprintf('%s\n%s', titlestr, descrip);
    end
    [hI, hCB] = dg_showGram(hF, f, f, corr, ...
        titlestr, ...
        '', '', 'covariance');
    dg_recolorGram(hCB, [-1 1], hI);
end

