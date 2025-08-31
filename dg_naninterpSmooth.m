function result = dg_naninterpSmooth(x, sw)
%result = dg_naninterpSmooth(x, sw)
% Linearly interpolates each run of NaNs in <x> with non-NaN values based
% on the non-NaN values surrounding the NaN run, smooths, and then sets the
% smoothed interpolated values back to NaN and returns the result.  If
% there is a NaN run at the beginning of <x>, the first non-NaN value is
% propagated backwards since interpolationg is impossible.  Similarly for
% the case where <x> ends with a NaN run.  Note that if there are non-NaN
% values at the endpoints, these do NOT propagate in the smoothing
% convolution which uses the Matlab default of padding with zeros on the
% ends.  If <x> is all NaN, then it is simply returned.  The gain is
% normalized so that a constant value remains the same (i.e.  the DC gain
% is unity).
%INPUTS
% x: the data series to be smoothed
% sw: the smoothing kernel.  I recommend the hanning(2*smoothing+1), where
%   <smoothing> is the number of points on one side of the center line.  In
%   any case, length(sw) must be an odd number to make the trimming of the
%   convolution work correctly.
%OUTPUTS
% result: the 'same' (i.e. central part of the) convolution; the 'valid'
%   part is given by result(smoothing + 1 : end - smoothing) where
%   <smoothing> is defined as "(length(sw) - 1) / 2".

%$Rev: 173 $
%$Date: 2013-04-24 16:05:36 -0400 (Wed, 24 Apr 2013) $
%$Author: dgibson $

if ~mod(length(sw), 2)
    error('dg_naninterpSmooth:swlength', ...
        '<sw> must contain an odd number of points');
end
swgain = sum(sw);
smoothing = (length(sw) - 1) / 2;
[nanidx nonnanidx] = dg_findNaNruns(x);
if isempty(nonnanidx)
    result = x;
    return
end
if isnan(x(1))
    % Propagate values into head NaNs & remove from <nanidx>
    x(1:nonnanidx(1)) = x(nonnanidx(1));
    headnans = 1:nonnanidx(1)-1;
    nanidx(1) = [];
else
    headnans = [];
end
if isnan(x(end))
    % Propagate values into tail NaNs & remove from <nanidx>
    x(nanidx(end):end) = x(nanidx(end)-1);
    tailnans = nanidx(end):length(x);
    nanidx(end) = [];
else
    tailnans = [];
end
% We now have nanidx and nonnanidx s.t. either nanidx is empty, or
% nonnanidx(1) < nanidx(1) and nonnanidx(end) > nanidx(end).  To make
% matched pairs of nanidx and nonnanidx s.t. nanidx(k) < nonnanidx(k), we
% remove the first nonnanidx:
nonnanidx(1) = [];
for runidx = 1:length(nanidx)
    x(nanidx(runidx):nonnanidx(runidx)-1) = interp1( ...
        [nanidx(runidx)-1, nonnanidx(runidx)], ...
        x([nanidx(runidx)-1, nonnanidx(runidx)]), ...
        nanidx(runidx):nonnanidx(runidx)-1);
end
s1 = conv(double(x), sw);
result = s1(smoothing + 1 : end - smoothing) / swgain;
for runidx = 1:length(nanidx)
    result(nanidx(runidx):nonnanidx(runidx)-1) = NaN;
end
result(headnans) = NaN;
result(tailnans) = NaN;

end

