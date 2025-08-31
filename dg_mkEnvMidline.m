function [midline, upper, lower] = dg_mkEnvMidline(waveform)
%INPUT
% waveform: a vector time series.
%OUTPUTS
% midline: calculated midline column vector.
% upper: upper envelope column vector.
% lower: lower envelope column vector.
%NOTES
% Envelopes are based on simple 3-point strict extrema, so they could miss
% extrema in integer-valued data where the top or bottom is strictly flat
% (strictly equal).

%$Rev: 311 $
%$Date: 2025-02-07 18:14:09 -0500 (Fri, 07 Feb 2025) $
%$Author: dgibson $

waveform = reshape(waveform, [], 1);
ispeak = [true 
    waveform(2:end-1)>waveform(1:end-2) & waveform(2:end-1)>waveform(3:end)
    true];
isvalley = [true 
    waveform(2:end-1)<waveform(1:end-2) & waveform(2:end-1)<waveform(3:end)
    true];
upper = interp1(find(ispeak), waveform(ispeak), (1:length(waveform))', ...
    'pchip');
lower = interp1(find(isvalley), waveform(isvalley), (1:length(waveform))', ...
    'pchip');
midline = mean([upper lower], 2);
