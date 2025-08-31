function f = dg_unwrap(x)
%f = dg_unwrap(x)  Unwraps a triangle wave of phase in units of radians.
%INPUT
% x: a vector of phase angles in radians, with values ranging from -pi/2 to
%   pi/2.  Behavior is undefined if any value exceeds that range.
%OUTPUT
% f: same as abs(diff(x)) except at extrema.  At the extrema, the
%   difference is calculated as if x were reflected about the closer of the
%   two values pi/2 or -pi/2.  "At" the extrema in this comment actually
%   means at the point following the extremum.

%$Rev: 53 $
%$Date: 2010-04-30 19:04:47 -0400 (Fri, 30 Apr 2010) $
%$Author: dgibson $

f = abs(diff(x));
peaks = findpeaks(x);
pts2fix = peaks + 1;
if ~isempty(peaks)
    if pts2fix(end) == length(x) + 1
        pts2fix(end) = [];
        peaks(end) = [];
    end
    f(pts2fix-1) = pi - x(pts2fix) - x(peaks);
end
troughs = findpeaks(-x);
if ~isempty(troughs)
    pts2fix = troughs + 1;
    if pts2fix(end) == length(x) + 1
        pts2fix(end) = [];
        troughs(end) = [];
    end
    f(pts2fix-1) = pi + x(pts2fix) + x(troughs);
end
    

% Copied from dg_emd, wherein they came from:
% http://www.mathworks.com/matlabcentral/fileexchange/19681-hilbert-huang-transform

function n = findpeaks(x)
% Find peaks.
% n = findpeaks(x)

n    = find(diff(diff(x) > 0) < 0);
u    = find(x(n+1) > x(n));
n(u) = n(u)+1;