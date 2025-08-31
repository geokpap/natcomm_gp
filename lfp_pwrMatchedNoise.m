function [result, units] = lfp_pwrMatchedNoise(fn, varargin)
%[result, units] = lfp_matchPowerNoise(fn)
% Function intended for use with lfp_createWave.
% Calculates FFT of the entire contents of filenum <fn>, then constructs a
% spectrum with magnitudes (and therefore also power) that match the
% original FFT but has uniformly randomly distributed phase.
%INPUTS
% fn: filenum
%OUTPUTS
% result: sample data of same size as lfp_Samples{fn}.
% units: same as lfp_SamplesUnits{fn}.
%OPTIONS
% 'enabledonly' - uses only data from enabled trials (see
%   lfp_enabledtrials).  Default is to use lfp_Samples{fn} in its entirety.
%   Since it is not easy to produce output containing a different number
%   of samples from the input, the noise samples are assigned to the same
%   samples in the output that were used as input, and the rest are set to
%   NaN.  This may cause trouble e.g. with filtering, in which case you
%   will have to replace the NaNs with something better.  Also note that
%   there will be discontinuities introduced at every point in the input
%   where there are non-enabled trials, which might affect the spectrum of
%   the output.
%NOTES
% If the signal you actually want to match contains big artifacts outside
% your analysis windows, you need to replace those big artifacts with
% signals that have the same spectrum as the good parts in order to use
% this function.
%   Regarding randomness: nothing is done here to initialize the random
% number generator, so it can be initialized as desired before calling this
% function.  Note that recent releases of Matlab (e.g. R2013a) give the
% same initial state to the random number generator when they are first
% started, which means that an uninitialized random number generator will
% generate the same sequence every time when using parProcess.
% Empirically, it seems "rng(0)" is the same as "rng('default')".

%$Rev: 349 $
%$Date: 2015-05-03 15:32:28 -0400 (Sun, 03 May 2015) $
%$Author: dgibson $

global lfp_Samples lfp_SamplesUnits lfp_TrialIndex

argnum = 1;
enabledonlyflag = false;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'enabledonly'
            enabledonlyflag = true;
        otherwise
            error('lfp_matchPowerNoise:badoption', ...
                ['The option "' varargin{argnum} '" is not recognized.'] );
    end
    argnum = argnum + 1;
end

if enabledonlyflag
    goodtrials = lfp_enabledTrials;
    isgoodsample = false(size(lfp_Samples{fn}));
    for tidx = 1:length(goodtrials)
        isgoodsample( lfp_TrialIndex(goodtrials(tidx), 3) : ...
            lfp_TrialIndex(goodtrials(tidx), 4) ) = true;
    end
    % We assume below that we have an even number of samples, so if
    % sum(isgoodsample) is odd, we throw out the very last one:
    if mod(sum(isgoodsample(:)), 2)
        lastone = find(isgoodsample, 1, 'last');
        isgoodsample(lastone) = false;
    end
else
    isgoodsample = true(size(lfp_Samples{fn}));
end

Y = reshape(fft(lfp_Samples{fn}(isgoodsample)), [], 1);
A = abs(Y);

% Construct noise spectrum by randomizing phase and enforcing symmetry
% rules for real signals.  If <A> is of even length, we avoid doing
% anything to the value exactly at Nyquist.  We rely here on the fact that
% the number of samples in lfp_Samples is always even, since CSCs are
% recorded in power-of-two-sized frames.
Nyqidx = round(length(A)/2) + 1;
A(2:Nyqidx-1) = A(2:Nyqidx-1) .* exp(1i*2*pi*rand(Nyqidx-2, 1));
A(Nyqidx+1:end) = conj(A((Nyqidx-1):-1:2));
result = NaN(size(lfp_Samples{fn}));
result(isgoodsample) = ifft(A);
if any(imag(result))
    error('lfp_matchPowerNoise:argh', ...
        'Matlab is broken!  Hardware failure!  Everyone to get from street!');
end
units = lfp_SamplesUnits{fn};

