function str = dg_dec2bin2(n)
%str = dg_dec2bin2(n)  Wrapper for dec2bin to convert negative binary
%   numbers as well as positive, by prepending a minus sign to the
%   converted value of -n.

%$Rev: 24 $
%$Date: 2009-03-31 21:51:08 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

if n < 0
    str = ['-' dec2bin(-n)];
else
    str = dec2bin(n);
end
