function str = dg_dec2bin(n,d)
%str = dg_dec2bin(n)
%   Returns a binary string representation of the floating point
%   number n, which must be between 0 and 1 inclusive.  <d> is optional and
%   specifies the number of digits (default 53, which is the precision of
%   a double float). 

%$Rev: 24 $
%$Date: 2009-03-31 21:51:08 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

if ~(n<=1 && n>=0)
    error('n must be between 0 and 1');
end

if nargin < 2
    d = 53;
end

nchars = d + 1;

str = repmat('.', 1, nchars);
rem = n;
exp = 1;
while (rem > 0 && exp < nchars)
    if rem >= 2^(-exp)
        rem = rem - 2^(-exp);
        str(exp+1) = '1';
    else
        str(exp+1) = '0';
    end
    exp = exp + 1;
end
for k = (exp+1):nchars
    str(k) = '0';
end
