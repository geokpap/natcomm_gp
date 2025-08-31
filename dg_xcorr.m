function [c,lags] = dg_xcorr(x,y,maxlag)
%dg_xcorr(x,y)
% Same as Matlab xcorr for real vector inputs, except that the scaling is
% done individually at each lag so that signals that are identical after
% performing the lag return a value of 1 at that lag.

%$Rev: 150 $
%$Date: 2012-05-09 17:52:40 -0400 (Wed, 09 May 2012) $
%$Author: dgibson $

if nargin < 3
    maxlag = max(length(x), length(y)) - 1;
end

if length(x) < length(y)
    x(end+1:length(y)) = 0;
elseif length(y) < length(x)
    y(end+1:length(x)) = 0;
end

lags = -maxlag:maxlag;
c = NaN(size(lags));
    
% Do the negative lags (i.e. shift y to the left):
for k = 1:maxlag
    c(k) = sum(x(1:end+lags(k)) .* y(-lags(k)+1:end)) / ... 
        sqrt(sum(x(1:end+lags(k)).^2) * sum(y(-lags(k)+1:end).^2));
end

% Do the non-negative lags (i.e. shift x to the left):
for k = maxlag+1:length(lags)
    c(k) = sum(x(lags(k)+1:end) .* y(1:end-lags(k))) / ... 
        sqrt(sum(x(lags(k)+1:end).^2) * sum(y(1:end-lags(k)).^2));
end
    
