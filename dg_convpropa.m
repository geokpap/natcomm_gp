function y = dg_convpropa(x, krnl)
% Your choice as to whether this function's name stands for "convolve with
% endpoint propagation" or "convolve propaly" (that's Bostonese).  It also
% normalizes by sum(krnl) so that the gain is unity at DC.  In other
% respects, it behaves like the 'same' option of 'conv'.
%INPUTS
% x: the series to be filtered as a vector of either orientation.
% krnl: the filter (smoothing) kernel as a column vector.
%OUTPUT
% y: the result of propagating the end points of <x> by as many points as
%   are needed to perform the convolution and return a vector the same
%   length as <x>.  I.e. the padding amounts to length(krn) - 1 points
%   at each end.
%NOTES
% The firs arg to 'conv' has length:
%   length(x) + 2*(length(krnl) - 1)
% so the length returned when using 'valid' is
%   length(x) + 2*(length(krnl) - 1) - length(krnl) + 1
%   = length(x) + 2*length(krnl) - 2 - length(krnl) + 1
%   = length(x) + length(krnl) - 1
% and to make it behave like the 'same' option means that we have to remove
% (length(krnl) - 1)/2 points at both ends, or if length(krnl) is even, we
% remove length(krnl)/2 from the left and length(krnl)/2 - 1 from the
% right.

%$Rev: 240 $
%$Date: 2016-03-18 16:27:51 -0400 (Fri, 18 Mar 2016) $
%$Author: dgibson $

pad = length(krnl) - 1;
y = conv( ...
    [x(1) * ones(pad, 1); reshape(x, [], 1); x(end) * ones(pad, 1)], ...
    krnl, 'valid' ) / sum(krnl);
y([(1:ceil((length(krnl) - 1)/2)) ...
    (end - floor((length(krnl) - 1)/2) + 1):end]) = [];

