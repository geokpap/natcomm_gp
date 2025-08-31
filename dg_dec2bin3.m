function bits = dg_dec2bin3(n)
%INPUT
% n: an integer to convert to a binary number represented as a bit
%   (Boolean) vector.
%OUTPUT
% bits: Boolean vector representation of <n> with bit 0 as the first
%   element.  Trailing zeros are suppressed.  Maximum length is 64.

%$Rev: 272 $
%$Date: 2021-02-12 12:46:53 -0500 (Fri, 12 Feb 2021) $
%$Author: dgibson $

bits = false(64,1);
rem = n;
if n == 0
    bits = false;
    return
end
while rem > 0
    bitidx = floor(log2(rem)) + 1;
    bits(bitidx) = true;
    rem = rem - 2^(bitidx - 1);
end
last1 = find(bits, 1, 'last');
bits(last1 + 1 : end) = [];

