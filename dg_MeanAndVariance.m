function result = dg_MeanAndVariance(A)

%dg_MeanAndVariance computes mean and variance from an array that may
%contain NaN values.
% result = dg_MeanRateAndVariance(A)
% result is an array with 2 rows and the same number of columns as rates;
% the first row contains the mean of each column vector in A after removing
% any elements that are NaN; the second row contains the "unbiased
% estimate" of the variance (i.e. using (N-1) in the denominator).  If a
% column of A contains only NaN values, then both rows in the same column
% of the return result are also NaN.

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

result = [];
for j = (1:length(A(1,:)))
    data = A(:,j);
    goodvalues = find(~isnan(data));
    if length(goodvalues) > 0
        newresults = [ mean(data(goodvalues)); var(data(goodvalues)) ];
    else
        newresults = [ NaN; NaN ];
    end
    result = [ result newresults ];
end