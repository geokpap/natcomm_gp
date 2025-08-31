function [binedges, binctrs] = dg_makeTimeBins(window, binwidth)
%[binedges, binctrs] = dg_makeTimeBins(window, binwidth, offsetflag)
% Computes integral-millisecond-width time bins that cover a window, using
% integer arithmetic until the final conversion from milliseconds to
% seconds.
%INPUTS
% window: time window in seconds to cover, where window(1) is the start
%   time and window(2) is the end time.  The actual time period covered
%   will be the closest that can be fit by an integral number of bins, so
%   the actual end of the window thus might be up to half a binwidth off in
%   either direction.  The start time of the window is exactly identical to
%   window(1).
% binwidth: in milliseconds.
%OUTPUTS
% binedges: in seconds
% binctrs: bin centers in seconds

%$Rev: 180 $
%$Date: 2013-09-10 17:31:17 -0400 (Tue, 10 Sep 2013) $
%$Author: dgibson $

if numel(window) > 2
    warning('dg_makeTimeBins:window', ...
        '<window> contains more than two elements.');
elseif numel(window) < 2
    error('dg_makeTimeBins:window2', ...
        '<window> must contain two elements.');
end
dur = round(diff(window) * 1000);
numbins = round(dur/binwidth);
binedges = (0 : numbins) * binwidth;
binctrs = (binedges(1:end-1) + binwidth/2) / 1000 + window(1);
binedges = binedges / 1000 + window(1);
