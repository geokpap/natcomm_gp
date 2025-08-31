function [f, P] = dg_sampleSpectrum(spectrum, freq, bw, verboseflag)
% Smooths spectrum with smoothing width proportional to center frequency,
% then finds frequency points with spacing inversely proportional to the
% magnitude of the second derivative of the smoothed spectrum.  Smoothing
% width is one point (i.e. no smoothing) at the first two freq points.  To
% keep it symmetrical, smoothing width is always an odd number of points.
% At the high frequency end, the last <halfwidth> points are linearly
% extrapolated from the preceding two.
%INPUTS
% spectrum: vector of log scaled spectral power.
% freq: vector of frequency points, one for each point in <spectrum>.  If
%   compatibility with dg_mknoise is desired, first frequency should be 0,
%   and the last frequency should be <Fs>/2.
% bw: the factor by which the frequency index is multiplied to calculate
%   the half-width of the smoothing kernel; must be less than 1.
% verboseflag: optional; should be <true> if given.
%OUTPUTS
% f: column vector of adaptively spaced frequency points suitable for input
%   to dg_mknoise.  Includes peak and valley points identified in the
%   smoothed version of <spectrum>.
% P: column vector of spectral power values transformed to linear scale
%   using P(f) = 10^(spectrum(f)/10).  With any luck, the result of
%       spline(f, sqrt(P), freq)
%   will be free of substantial overshoots of the original <spectrum>.

%$Rev: 266 $
%$Date: 2019-08-28 17:07:36 -0400 (Wed, 28 Aug 2019) $
%$Author: dgibson $

if nargin < 4
    verboseflag = false;
end
if numel(freq) ~= numel(spectrum)
    error('dg_sampleSpectrum:numel', ...
        '<spectrum> and <freq> must contain the same number of elements.');
end
if bw >= 1
    error('dg_sampleSpectrum:bw', ...
        '<bw> must be less than 1.');
end

smoothspec = zeros(numel(spectrum), 1);
prevhalfwidth = NaN;
for fidx = 1:length(freq)
    halfwidth = floor((fidx - 1) * bw);
    if halfwidth ~= prevhalfwidth
        if verboseflag
            fprintf('Starting halfwidth=%d\n', halfwidth);
        end
        hw = hanning(2 * halfwidth + 1);
    end
    if fidx + halfwidth <= length(freq)
        smoothspec(fidx) = 10^(( ...
            sum(hw .* spectrum(fidx + (-halfwidth:halfwidth))) / sum(hw) ...
            ) / 10);
    else
        smoothspec(fidx) = smoothspec(fidx-1) * ...
            smoothspec(fidx-1)/smoothspec(fidx-2);
    end
end

maxidx = find( smoothspec(2:end-1) >= smoothspec(1:end-2) & ...
    smoothspec(3:end) < smoothspec(2:end-1) );
minidx = find( smoothspec(2:end-1) < smoothspec(1:end-2) & ...
    smoothspec(3:end) >= smoothspec(2:end-1) );
fidx = sort([maxidx; minidx]);
% Add another point to the end at equal spacing to the last two:
fidx(end+1) = fidx(end) + diff(fidx(end-1:end));
% Include the first and last freq points:
fidx = [1; fidx; length(freq)];
for k = 1:2
    fidx = unique([fidx; fidx(1:end-1) + round(diff(fidx)/2)]);
end
% <intrpthresh> is the ratio of the difference-of-slopes to the smoothed
% spectrum that triggers interpolation of another frequency point:
intrpthresh = 1; 
needsintrp = 0;
while ~isempty(needsintrp)
    % number of points between each pair in <fidx>:
    numpts = diff(fidx);
    % slope is amplitude change per <freq> point:
    slope = diff(sqrt(smoothspec(fidx))) ./ numpts;
    dslope = diff(slope);
    % dslope(n) is the difference in slope around fidx(n+1). When dslope(n)
    % is large, that means that points fidx(n+1) and fidx(n+2) need
    % interpolation.  "Large" means relative to the power values, but
    % weighted by the number of points (<numpts>) in the next piecewise
    % linear bit of the approximated spectrum. <needsintrp> is an index
    % into fidx for the center point that needs interpolation on the right.
    needsintrp = find(abs(numpts(2:end) .* dslope ./ sqrt(smoothspec(fidx(3:end)))) ...
        > intrpthresh ) + 1;
    % Interpolate points needsintrp and needsintrp+1. This formula
    % is simplified from: sqrt( ...
    %   fidx(needsintrp+1) ./ fidx(needsintrp) ) .* fidx(needsintrp)
    newfpts = round(sqrt( ...
        fidx(needsintrp+1) .* fidx(needsintrp) ));
    if all(ismember(newfpts, fidx))
        % no more interpolations possible! (or <needsintrp> is empty)
        needsintrp = [];
    else
        fidx = unique([fidx; newfpts]);
    end
end
f = freq(fidx);
P = smoothspec(fidx);
