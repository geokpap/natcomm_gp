function [S, freqedges, f, A] = dg_hht(imf, varargin)
%S = dg_hht(imf, h, f) Empirical Mode Decomposition (Hilbert-Huang Transform)
% Produces a matrix suitable for display as a spectrogram.
% NaN values in f propagate through fsmoothing and into the final gram.
%INPUTS
% imf: an array containing IMFs in components x points form.  The "residual
%   junk" component that remains after extracting the true IMFs should not
%   be included in <imf>.
%OUTPUTS
% S: the matrix that corresponds to a spectrogram, with one column for each
%   time point in <imf> and one row for each frequency bin.  Each time-freq
%   point represents the average power over all IMFs.
% freqedges: the edges that define the frequency bins for <S>.  Note that
%   these will be logarithmically spaced if 'logf' is specified.
% f: instantaneous frequency traces, same size as imf
% A: as returned by dg_AMFMdecomp.
%OPTIONS
% 'dq' - use Direct Quadrature method to calculate instantaneous
%   frequencies from normalized IMFs returned by dg_AMFMdecomp (default
%   method is Normalized Hilbert Transform) 
% 'f', f - specify the instantaneous frequencies that go with the IMFs
%   instead of calculating them from <imf>.  As returned by dg_emd, in
%   components x points form.  See 'h' option.
% 'fblur', n - blur the final gram by <n> points in the frequency
%   direction. Blurring is done with a Hanning window 2*n+1 points wide.
% 'freqedges', freqedges - specify the edges of the frequency bins in
%   radians per sample; default is <numfreqs> bins spanning the range of
%   values in the smoothed frequency traces, where <numfreqs> is specified
%   by 'numfreqs' option.
% 'fsmoothing', n - smooth the frequency traces calculated from the EMD
%   components with a Hanning window 2*n+1 points wide.  This affects both
%   of the return values <S> and <f>, as well as the plots produced by
%   'plotfreqs'.
% 'h', h - specify the Hilbert Analytic Signals for the IMFs instead of
%   calculating them from <imf>.  As returned by dg_emd, in components x
%   points form.  If 'h' is specified without specifying 'f', then the
%   frequencies are calculated by differentiating the phases of <h>.
% 'logf' - bin frequencies on a log scale (default is linear)
% 'method', method - for dg_AMFMdecomp, ultimately passed through to Matlab
%   function 'interp1'.  Default: 'spline'.
% 'numfreqs', n - use <n> bins in the frequency dimension; default = 1000.
% 'plotfreqs', figtitle, Ts, startpt - creates a figure with a plot of
%   instantaneous freq vs. time for all IMFs, titled <figtitle>;
%   axes are calibrated based on the sample period <Ts> given in seconds.
%   <startpt> is the absolute point number of the first time point, i.e.
%   times are displayed as Ts * ((1:size(imf,2) - 1) - startpt.
% 'tblur', n - like 'fblur', but operates in the time direction.
% 'trimming', trimming -  Works like the <trimming> parameter to the
%   dg_emd 'hilbert' option, except that it looks at the first component
%   instead of the raw waveform.  If <trimming> is 'extrema', then the
%   frequencies in the interval from the beginning through the first
%   maximum or minimum are set to NaN, and likewise from the last maximum
%   or minimum to the end. If trimming is a number, then this number of
%   samples is set to NaN at each end.  If it is a two-element vector, then
%   the first element is the number of NaNs at the beginning and the second
%   is the number at the end.  Default is to do no trimming, but the entire
%   first column of S is NaN anyway because instantaneous freq (being the
%   first difference of phase) is not defined at the first sample.
% 'whiten' - multiply <S> by nominal bin frequencies to whiten pink
%   signals; the nominal frequencies are the averages of the <freqedges>
%   that bound each bin.

%$Rev: 60 $
%$Date: 2010-07-20 17:42:35 -0400 (Tue, 20 Jul 2010) $
%$Author: dgibson $

A=[];

dqflag = false;
f = [];
fblur = 0;
figtitle = NaN;
freqedges = [];
fsmoothing = 0;
h = [];
logfflag = false;
method = 'spline';
numfreqs = 1000;
tblur = 0;
trimming = [];
whitenflag = false;
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'dq'
            dqflag = true;
        case 'f'
            argnum = argnum + 1;
            f = varargin{argnum};
        case 'fblur'
            argnum = argnum + 1;
            fblur = varargin{argnum};
            if ~isnumeric(fblur)
                error('dg_hht:badfblur', ...
                    '''fblur'' requires a numeric argument');
            end
        case 'freqedges'
            argnum = argnum + 1;
            freqedges = varargin{argnum};
            if ~isnumeric(fsmoothing)
                error('dg_hht:badfreqedges', ...
                    '''freqedges'' requires a numeric argument');
            end
            numfreqs = length(freqedges) - 1;
        case 'fsmoothing'
            argnum = argnum + 1;
            fsmoothing = varargin{argnum};
            if ~isnumeric(fsmoothing)
                error('dg_hht:badfsmoothing', ...
                    '''fsmoothing'' requires a numeric argument');
            end
        case 'h'
            argnum = argnum + 1;
            h = varargin{argnum};
        case 'logf'
            logfflag = true;
        case 'method'
            argnum = argnum + 1;
            method = varargin{argnum};
        case 'numfreqs'
            argnum = argnum + 1;
            if isempty(freqedges)
                numfreqs = varargin{argnum};
            end
        case 'plotfreqs'
            argnum = argnum + 1;
            figtitle = varargin{argnum};
            argnum = argnum + 1;
            Ts = varargin{argnum};
            argnum = argnum + 1;
            startpt = varargin{argnum};
            if ~isnumeric(Ts) || ~isscalar(Ts)
                error('dg_hht:Ts', ...
                    'The <Ts> arg for ''plotfreqs'' must be a number');
            end
            if ~ischar(figtitle)
                error('dg_hht:figtitle', ...
                    'The <figtitle> arg for ''plotfreqs'' must be a string');
            end
            if ~isnumeric(startpt) || ~isscalar(startpt)
                error('dg_hht:startpt', ...
                    'The <Ts> arg for ''startpt'' must be a number');
            end
        case 'tblur'
            argnum = argnum + 1;
            tblur = varargin{argnum};
            if ~isnumeric(tblur)
                error('dg_hht:badtblur', ...
                    '''tblur'' requires a numeric argument');
            end
        case 'trimming'
            argnum = argnum + 1;
            trimming = varargin{argnum};
            if ~( isequal(trimming, 'extrema') || isnumeric(trimming) )
                error('dg_hht:trimming', ...
                    '"%s" is not a legal value for <trimming>', ...
                    dg_thing2str(trimming) );
            end
        case 'whiten'
            whitenflag = true;
        otherwise
            error('dg_hht:badoption', ...
                ['The option "' dg_thing2str(varargin{argnum}) '" is not recognized.'] );
    end
    argnum = argnum + 1;
end

if isempty(f)
    if isempty(h)
        if dqflag
            for cidx = 1:size(imf,1)
                [A(cidx,:), F(cidx,:)] = ...
                    dg_AMFMdecomp(imf(cidx,:), method);
                th = atan(F(cidx,:) ./ sqrt(1 - F(cidx,:).^2)); % radians
                f(cidx,:) = [NaN diff(dg_unwrap(th))]; % radians per sample
            end
        else
            for cidx = 1:size(imf,1)
                [f(cidx,:), A(cidx,:), F(cidx,:)] = ...
                    dg_nht(imf(cidx,:), method);
            end
        end
        P = A.^2;   % power
    else
        th = unwrap(angle(h')); % radians
        f = [NaN(size(th,2)); diff(th)]'; % radians per sample
        P = conj(h) .* h;   % power
        A = abs(h);
    end
end

if any(f(:)<0)
    msg = sprintf('f<0 at %d samples', sum(f(:)<0) );
    for cidx = 1:size(imf,1)
        if any(f(cidx,:)<0)
            msg = sprintf( ...
                '%s\nimf %d: %d samples', msg, cidx, sum(f(cidx,:)<0) );
        end
    end
    warning('dg_hht:negf', '%s', msg');
end
if any(A(:)<0)
    msg = sprintf('A<0 at %d samples', sum(A(:)<0) );
    for cidx = 1:size(imf,1)
        if any(A(cidx,:)<0)
            msg = sprintf( ...
                '%s\nimf %d: %d samples', msg, cidx, sum(A(cidx,:)<0) );
        end
    end
    warning('dg_hht:negA', '%s', msg');
end

if ~isempty(trimming)
    if isequal(trimming, 'extrema')
            peaks = dg_findpks(imf(1,:));
            troughs = dg_findpks(-imf(1,:));
            if isempty(peaks) || isempty(troughs) ...
                    || peaks(1) == 1 || troughs(1) == 1 ...
                    || peaks(end) == size(imf,2) || troughs(end) == size(imf,2)
                error('dg_hht:oops', ...
                    'Pathological case, panic!');
            end
            startidx = max(peaks(1), troughs(1)) + 1;
            endidx = min(peaks(end), troughs(end)) - 1;
    else
        startidx = trimming(1);
        if isscalar(trimming)
            endidx = trimming;
        else
            endidx = trimming(2);
        end
    end
    f(:, [1:startidx-1 endidx+1:end]) = NaN;
    P(:, [1:startidx-1 endidx+1:end]) = NaN;
end


% smooth the freq traces, excluding the terminal NaNs
hw = hanning(2*fsmoothing+1);
for c = 1:size(imf,1)
    non_nans = ~isnan(f(c,:));
    smoothed = conv(f(c,non_nans), hw) / sum(hw);
    f(c,non_nans) = smoothed(fsmoothing + 1 : end - fsmoothing);
end

if ~isequalwithequalnans(figtitle, NaN)
    dg_hht_plotfreqs(f, A, Ts, startpt, figtitle);
end

% Compute time-freq matrix with same number of cols as x has points,
% containing the Hilbert analytic function power averaged over all IMFs.
% Like histc, uses freqedges(k) <= f < freqedges(k+1).
% <binned_P> is in freq x sample x IMF format.
binned_P = zeros(numfreqs, size(imf,2), size(imf,1));
maxfreq = max(f(:));
if isempty(freqedges)
    if logfflag
        % Handle non-positive freq points:
        if any(any(f<=0))
            f(f<=0) = min(f(f>0));
            warning('dg_hht:badfreqs', ...
                'Replaced non-positive smoothed freqs with smallest positive value');
        end
        minfreq = min(f(:));
        freqedges = logspace( floor(log10(minfreq)), ceil(log10(maxfreq)), ...
            numfreqs+1 );
    else
        minfreq = 0;
        freqedges = linspace(minfreq, maxfreq, numfreqs+1);
    end
end
for imfidx = 1:(size(imf,1))
    for k = 1:numfreqs
        thepoints = freqedges(k) <= f(imfidx,:) ...
            & f(imfidx,:) < freqedges(k+1);
        binned_P(k, thepoints, imfidx) = P(imfidx, thepoints);
    end
end
if fblur || tblur
    hw2 = repmat(hanning(2*fblur+1), 1, 2*tblur+1) .* ...
        repmat(hanning(2*tblur+1)', 2*fblur+1, 1);
    S = conv2(squeeze(mean(binned_P, 3)), hw2);
    S = S(fblur+1:end-fblur, tblur+1:end-tblur);
else
    S = squeeze(mean(binned_P, 3));
end

% Make NaNs propagate
hasnan = any(isnan(f));
S(:, hasnan) = NaN;

if whitenflag
    freqs = (freqedges(1:end-1) + (freqedges(2:end))/2);
    S = S .* repmat(reshape(freqs, [], 1), 1, size(S, 2));
end


