function [M, N] = dg_nanTolerantMedian(A, dim)
%[M, N] = dg_nanTolerantMedian(A, dim)
% Works for 2-D arrays exactly like Matlab 'sum', except that NaNs are
% ignored.  The mean of a vector of all NaNs is defined to be NaN.  N is
% the number of non-NaNs that contributed to each value in M.
%
% NOTES: 
% 1) dim is optional, and if given, it must be 1 or 2.
% 2) Use Matlab nanmedian if N is not required.

%$Rev: 24 $
%$Date: 2009-03-31 21:51:08 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

if nargin < 2
    dim = 1;
end

if dim == 2
    A = A';
    M = NaN(size(A,2), 1);
else
    M = NaN(1, size(A,2));
end

N = M;
for k = 1:size(A,2)
    N(k) = sum(~isnan(A(:,k)));
    if N(k)
        M(k) = median(A(~isnan(A(:,k)), k), 1);
    end
end

