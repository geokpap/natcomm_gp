function [result, units] = lfp_matchPowerNoise(fn)
%[result, units] = lfp_matchPowerNoise(fn)
% Function intended for use with lfp_createWave.
% Calculates FFT of the entire contents of filenum <fn>, then constructs a
% spectrum with magnitudes (and therefore also power) that match the
% original FFT but has uniformly randomly distributed phase.
%INPUTS
% fn: filenum
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

%$Rev: 343 $
%$Date: 2015-04-05 16:41:55 -0400 (Sun, 05 Apr 2015) $
%$Author: dgibson $

global lfp_Samples lfp_SamplesUnits

Y = reshape(fft(lfp_Samples{fn}(:)), [], 1);
A = abs(Y);

% Construct noise spectrum by randomizing phase and enforcing symmetry
% rules for real signals:
Nyqidx = length(A)/2 + 1;
A(2:Nyqidx-1) = A(2:Nyqidx-1) .* exp(1i*2*pi*rand(Nyqidx-2, 1));
A(Nyqidx+1:end) = conj(A((Nyqidx-1):-1:2));
result = ifft(A);
if any(imag(result))
    error('lfp_matchPowerNoise:argh', ...
        'Matlab is broken!  Hardware failure!  Everyone to get from street!');
end
units = lfp_SamplesUnits{fn};

