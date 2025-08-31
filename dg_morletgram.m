function [P, freqs, wt] = dg_morletgram( ...
    x, dt, fspacing, freqlim, param, fscale)
%INPUTS
% x: time series vector.
% dt: sample period.
% fspacing: ratio between scales on successive rows.
% freqlim: two-element vector where freqlim(1) is the the pseudo-frequency
%   of the longest scale and freqlim(2) is a limit which will not be
%   exceeded by the pseudo-frequency of the shortest scale.
% param: optional argument to specify the Morlet wavelet width parameter.
%   This is approximately the number of cycles in the wavelet. Default = 6.
%   Values below about 4 or 5 are considered dangerous, whereas above about
%   15 they become boring.
% fscale: string to specify whether the pseudo-frequencies should be spaced
%   linearly ('lin') or logarithmically ('log').  Default = 'log'.  When
%   using 'lin', <fspacing> specifies the difference between successive
%   pseudo-frequencies (scales) rather than there ratio.
%OUTPUTS
% P: pseudo-power, i.e. the squared magnitude of <wt.cfs>.
% freqs: vector of pseudo-frequency values corresponding to the scale on
%   each row of <wt>.
% wt: the complex continuous wavelet transform of <x>, as returned by
%   cwtft.

%$Rev: 222 $
%$Date: 2015-07-09 20:52:29 -0400 (Thu, 09 Jul 2015) $
%$Author: dgibson $

if nargin < 5
    param = 6;
end
if nargin < 6
    fscale = 'log';
end

switch fscale
    case 'log'
        % convert everything to log units to construct the frequency grid.
        logfspacing = log(fspacing);
        logfreqlim = log(freqlim);
        logf = logfreqlim(1) : logfspacing : logfreqlim(2);
        freqs = exp(logf);
    case 'lin'
        freqs = freqlim(1) : fspacing : freqlim(2);
    otherwise
        error('dg_morletgram:fscale', ...
            'Bad value for <fscale>: %s', fscale);
end

sig1 =  struct('val',x,'period',dt);
FourierFactor = 4*pi/(param+sqrt(2+param^2));
scales = 1./(FourierFactor * freqs);
wt = cwtft(sig1, 'scales', scales, 'wavelet', {'morl', param});
P = wt.cfs .* conj(wt.cfs);

