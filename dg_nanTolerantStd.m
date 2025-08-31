function [s, N] = dg_nanTolerantStd(A, flag, dim)
%[s, N] = dg_nanTolerantStd(A, dim)
% Works for 2-D arrays exactly like Matlab 'std', except that NaNs are
% ignored.  The mean of a vector of all NaNs is defined to be NaN.  N is
% the number of non-NaNs that contributed to each value in s.
%
% NOTES: 
% 1) dim is optional, and if given, it must be 1 or 2.
% 2) Use Matlab nanstd if N is not required.
% 3) flag is as for std.

%$Rev: 24 $
%$Date: 2009-03-31 21:51:08 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

if nargin < 2
    flag = 0;
end
if nargin < 3
    dim = 1;
end

if dim == 2
    A = A';
    s = NaN(size(A,2), 1);
else
    s = NaN(1, size(A,2));
end

N = s;
for k = 1:size(A,2)
    N(k) = sum(~isnan(A(:,k)));
    if N(k)
        s(k) = std(A(~isnan(A(:,k)), k), flag, 1);
    end
end

