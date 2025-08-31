function [Pxx, f, SEM, numwin] = dg_pwelch(x, winlen, noverlap, Fs, varargin)
% Similar to a subset of 'pwelch', but also offers winlen-by-winlen
% detrending, and scaling so that sum(Pxx) = 1.  The taper is always a
% hanning window.  The first window begins at the first point in <x>, and
% up to <winlen> - <noverlap> - 1 points may be ignored at the end of <x>.
%INPUTS
% x: time series to transform.
% winlen: number of time points in each data window.
% overlap: number of samples by which to overlap successive windows.  May
%   be [], in which case round(length(x)/2) is used.
% Fs: sampling rate.
%OUTPUT
% Pxx: average power spectrum.
% f: frequencies of each point in <P>.
% SEM: standard error of mean at each point in <P>.
% numwin: the number of windows that went into the calculations.
%OPTIONS
% 'detrend' - removes trend from each winlen before applying the taper.
% 'norm' - scales so that sum(Pxx) = 1.  <SEM> is scaled by the same
%   factor.
% 'pad', N - as for lfp_mtspectrum.  Unless window width happens already to
%   be 2^k samples long for some k, there is always padding, and N controls
%   how many powers of 2 to skip: N=0 pads to the next greater power of 2,
%   N=1 pads to twice that length, N=2 pads to four times that length, etc.
%   The default value is N=0.


%$Rev:  $
%$Date:  $
%$Author: dgibson $

detrendflag = false;
normflag = false;
pad = 0;

argnum = 0;
while true
    argnum = argnum + 1;
    if argnum > length(varargin)
        break
    end
    if ~ischar(varargin{argnum})
        continue
    end
    switch varargin{argnum}
        case 'detrend'
            detrendflag = true;
        case 'norm'
            normflag = true;
        case 'pad'
            argnum = argnum + 1;
            pad = varargin{argnum};
        otherwise
            error('dg_pwelch:badoption', ...
                'The option %s is not recognized.', ...
                dg_thing2str(varargin{argnum}));
    end
end
if isempty(noverlap)
    noverlap = round(length(x)/2);
end
if noverlap >= winlen
    error('dg_pwelch:noverlap', ...
        '<noverlap> must be less than <winlen>.');
end
exponent = nextpow2(winlen) + pad;
nfft = 2^exponent;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This bit is copied from Chronux.  I have faith.
df = Fs/nfft;
f = 0:df:Fs; % all possible frequencies
f = f(1:nfft);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
winshift = winlen - noverlap;
numwin = floor((length(x) - winlen) / winshift) + 1;
datary = NaN(winlen, numwin);
for winidx = 1:numwin
    datary(:, winidx) = x((1:winlen) + (winidx - 1) * winshift);
end
if detrendflag
    datary = detrend(datary);
end
taper = hanning(winlen);
S = fft(datary .* repmat(taper, 1, numwin), nfft);
P = conj(S) .* S;

BL.sum = sum(P, 2);
BL.sumsqrs = sum(P.^2, 2);
Pxx = BL.sum/numwin;
SEM = sqrt((BL.sumsqrs - (BL.sum).^2/numwin)/((numwin-1)*numwin));

if normflag
    SEM = SEM/sum(Pxx);
    Pxx = Pxx/sum(Pxx);
end
