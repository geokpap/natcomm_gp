function n = dg_findpks(x)
%n = dg_findpks(x)
% This is rewritten from the "findpeaks" function from
% http://www.mathworks.com/matlabcentral/fileexchange/19681-hilbert-huang-transform
% Its strength is that it simply defines peaks as points where x ceases to
% be monotonically increasing, so there is no issue about how to handle
% plateaus. Its weakness is the same thing.
%OUTPUT
% Returns a row vector.

%$Rev: 89 $
%$Date: 2010-11-01 18:06:05 -0400 (Mon, 01 Nov 2010) $
%$Author: dgibson $

isincreasing = diff(x(:)) > 0;
n = reshape(find(isincreasing(1:end-1) & ~isincreasing(2:end)) + 1, 1, []);
