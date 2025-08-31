function result = dg_nanTolerantSum(A)
%Deprecated.  Use Matlab nansum instead.
%dg_nanTolerantSum computes the sum of an array that may
%contain NaN values.
% result = dg_MeanRateAndVariance(A)
% For 2-D arrays, this works just like the Matlab sum function, except that
% NaN values are ignored in computing the sum.  If a column of A contains
% only NaN values, then the same column of the return result is also NaN.

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

result = [];
for j = (1:length(A(1,:)))
    data = A(:,j);
    goodvalues = find(~isnan(data));
    if length(goodvalues) > 0
        newresult = sum(data(goodvalues));
    else
        newresult = NaN;
    end
    result = [ result newresult ];
end