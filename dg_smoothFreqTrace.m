function f = dg_smoothFreqTrace(freq, smoothing)
%f = dg_smoothFreqTrace(freq, smoothing) does the exact same smoothing as
% is invoked by the 'fsmoothing' option in dg_nht.
%INPUTS
% freq: a series of time samples whose first element is worthless.
% smoothing: smooths with a Hanning window 2*<smoothing>+1
%   points wide, normalized for unity gain at DC.  The initial NaN is
%   removed before the smoothing and then replaced afterwards. Matlab's
%   'conv' implicitly extends f with zeros at the beginning and end, so
%   counting the initial NaN, there are <width>+1 samples at the beginning
%   and <width> at the end that are "invalid" (i.e. either NaN or affected
%   by the added zeros). 

%$Rev: 96 $
%$Date: 2011-01-06 21:15:40 -0500 (Thu, 06 Jan 2011) $
%$Author: dgibson $

hw = hanning(2*smoothing+1);
f = conv(hw, freq(2:end)) / sum(hw);
f = [NaN f(smoothing+1:end-smoothing)];
