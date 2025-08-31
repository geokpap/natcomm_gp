function [c,lags] = dg_xcov(x,y,maxlag)
%dg_xcorr(x,y)
% Same as Matlab xcov for real vector inputs, except that the scaling is
% done individually at each lag so that signals that are identical after
% performing the lag return a value of 1 at that lag.
%NOTES
% A crude benchmark consisting of running either dg_xcov(w(:,k), x(:,k), 0)
% or xcov(w(:,k), x(:,k), 0, 'coeff') on each column k of two 100 x 1000
% matrices showed dg_xcov to run about five times faster than xcov.
% However, setting nlags = 20 (i.e., running dg_xcov(w(:,k), x(:,k), 20))
% showed dg_xcov to be about 1.4 times slower than xcov.

%$Rev: 162 $
%$Date: 2012-12-15 20:48:35 -0500 (Sat, 15 Dec 2012) $
%$Author: dgibson $

if nargin < 3
    maxlag = max(length(x), length(y)) - 1;
end

[c,lags] = dg_xcorr(x-mean(x),y-mean(y),maxlag);

