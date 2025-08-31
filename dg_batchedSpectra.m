function [Spectra, f] = dg_batchedSpectra(data, varargin)
%[Spectra, f] = dg_batchedSpectra(data)
% Applies multi-taper spectral analysis to each column of <data>.
% Same as dg_crossFreqCorr Rev 36, except returns Spectra instead of doing
% cross-freq correlation.  Does not plot anything.
%INPUTS
%  <data> contains one time window of consecutive samples in each column.
%  The numbers of rows (samples per window) and columns (windows) are both
%  variable.
%OUTPUTS
%  <Spectra> contains the power spectrum for each window, in [freq X
%   window] format. 
%  <f> is the vector of frequencies represented by the matrix.
%OPTIONS
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
%  be analyzed in batches.  In this case, <data> does NOT contain the
%  actual data, but is a 2-column array that contains starting point
%  indices in the first column and ending point indices in the second
%  column, constrained s.t. the number of points between the starting point
%  and ending point is always the same.  The indices point into
%  <datasource>.  Neither <datasource> nor <data> are altered in any way,
%  so Matlab "should" pass them by reference (i.e., not copy them).
%'k', numtapers
%  Forces number of multitapers to <numtapers>; default is 1 (note that
%  this differs from lfp_mtspectrum).
%'nw', nw
%  Forces "time-bandwidth product" setting for multitapers to <nw>; default
%  is 1.8 (note that this differs from lfp_mtspectrum).
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

%$Rev: 197 $
%$Date: 2014-05-29 17:22:15 -0400 (Thu, 29 May 2014) $
%$Author: dgibson $

datasource = [];
descrip = '';
detrendflag = true;
freqlim = [];
fs = 0;
K = 1;
NW = 1.8;
pad = 0;
plotflag = true;
verboseflag = false;
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'detrend'
            argnum = argnum + 1;
            detrendflag = varargin{argnum};
            if ~detrendflag
                rmdcflag = false;
            end
        case 'freqlim'
            argnum = argnum + 1;
            freqlim = varargin{argnum};
        case 'fs'
            argnum = argnum + 1;
            fs = varargin{argnum};
        case 'getfrom'
            argnum = argnum + 1;
            datasource = varargin{argnum};
            if isempty(datasource)
                error('dg_crossFreqCorr:datasource', ...
                    '<datasource> is empty');
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
        case 'nw'
            argnum = argnum + 1;
            NW = varargin{argnum};
        case 'pad'
            argnum = argnum + 1;
            if ~strcmp(class(varargin{argnum}), 'double')
                error('lfp_spec:badpad', ...
                    'The padding factor must be a number');
            end
            pad = varargin{argnum};
        case 'rmdc'
            argnum = argnum + 1;
            rmdcflag = varargin{argnum};
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
            memory_array(i) = mem_elements(i).bytes;
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
    disp(sprintf('Processing in approx. %d batches', estbatches));
end

% First pass
done = false;
sumP = zeros(numfreqs2use, 1);
sumsqrsP = zeros(numfreqs2use, 1);
N = 0;
batchnum = 0;
tempfilenames = {};
while ~done
    batchnum = batchnum + 1;
    if verboseflag
        disp(sprintf( ...
            'Batch %.0f, windows %.0f-%.0f', ...
            batchnum, startwin, endwin ));
    end
    % compute one batch of power spectra
    numwin = endwin - startwin + 1;
    if isempty(datasource)
        batchsamples = data(:, startwin:endwin);
    else
        samplidx = repmat((0:winsamples-1)', 1, numwin);
        startidx = repmat(data(startwin:endwin,1)', winsamples, 1);
        batchsamples = datasource(startidx + samplidx);
    end
    clear samplidx startidx
    if isequal(detrendflag, true)
        batchsamples = detrend(batchsamples);
    elseif isequal(detrendflag, 'rmdc')
        batchsamples = batchsamples - ...
            repmat(mean(batchsamples,1), size(batchsamples,1), 1);
    end
    % S comes back in form frequency index x taper index x channels/trials:
    S = chronux_mtfftc(batchsamples, tapers, nfft);
    clear batchsamples
    % Compute <P>  (spectral power), in [freq X window] format - i.e. each
    % column of <P> represents the same column of <data>.
    P = squeeze(mean(conj(S(isfreq2use,:,:)) .* S(isfreq2use,:,:), 2));
    clear S
    sumP = sumP + sum(P,2);
    sumsqrsP = sumsqrsP + sum(P.^2,2);
    N = N + size(P, 2);
    if ismultibatch
        tempfilenames{end+1} = sprintf('%s%07d.mat', tempbase, startwin);
        if verboseflag
            disp(sprintf( ...
                'Saving file %s', ...
                tempfilenames{end} ));
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
avgP = sumP/numwins;
% We use the sqrt-sum-of-squared-diffs, i.e. sqrt(N-1)*std, instead of
% std for sigma, which must be what Masimore ... Redish 2004 eqn 1
% meant, because otherwise the final value is proportional to N-1.  Using
% varP = (sumsqrsP - (sumP).^2/N)/(N-1), we get the formula for sigma:
sig = sqrt(sumsqrsP - sumP.^2/numwins);
    
% Second pass
done = false;
xcovsum = zeros(numfreqs2use);
if ismultibatch
    numbatches = length(tempfilenames);
else
    numbatches = 1;
end
Spectra = zeros(size(P,1), 0);
for tempfilenum = 1:numbatches
    if ismultibatch
        if verboseflag
            disp(sprintf( ...
                'Loading temp file #%.0f: %s', ...
                tempfilenum, tempfilenames{tempfilenum} ));
        end
        load(tempfilenames{tempfilenum});
        delete(tempfilenames{tempfilenum});
    end
    Spectra = [ Spectra P ];
end    



