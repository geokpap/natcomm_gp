function [f, A, F, h] = dg_nht(imf, method, varargin)
%[f, A, F] = dg_nht(imf)  Normalized Hilbert Transform
%INPUTS
% imf: an individual IMF as A row vector.
% method: for dg_AMFMdecomp
%OUTPUTS
% f: instantaneous frequency in radians per sample between points k-1 and k
%   for k>2.  First element is NaN.  (Nyquist = pi rad/sample.)
% A: amplitude component from dg_AMFMdecomp.
% F: frequency component from dg_AMFMdecomp.
% h: hilbert(F).
%OPTIONS
% 'fsmoothing', width - smooths <f> with a Hanning window 2*<width>+1
%   points wide, normalized for unity gain at DC.  The initial NaN is
%   removed before the smoothing and then replaced afterwards. Matlab's
%   'conv' implicitly extends f with zeros at the beginning and end, so
%   counting the initial NaN, there are <width>+1 samples at the beginning
%   and <width> at the end that are "invalid" (i.e. either NaN or affected
%   by the added zeros). 

%$Rev: 95 $
%$Date: 2010-12-31 21:59:34 -0500 (Fri, 31 Dec 2010) $
%$Author: dgibson $

%REFERENCES
% Huang NE, Wu Z, Long SR, Arnold KC, Chen X, Blank K, "On Instaneous
% Frequency", Advances in Adaptive Data Analysis Vol. 1, No. 2 (2009)
% 177?229

fsmoothing = 0;
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'fsmoothing'
            argnum = argnum + 1;
            fsmoothing = varargin{argnum};
        otherwise
            error('funcname:badoption', ...
                ['The option "' dg_thing2str(varargin{argnum}) '" is not recognized.'] );
    end
    argnum = argnum + 1;
end

[A, F] = dg_AMFMdecomp(imf, method);
h = hilbert(F);
th = unwrap(angle(h));
if fsmoothing == 0
    f = [NaN diff(th)];
else
    hw = hanning(2*fsmoothing+1);
    f = conv(hw, diff(th)) / sum(hw);
    f = [NaN f(fsmoothing+1:end-fsmoothing)];
end
