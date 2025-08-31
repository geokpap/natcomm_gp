function [numglitch, hF] = dg_superglitchviewer(samples, thresh, maxpts)
%hF = dg_superglitchviewer(samples, thresh, maxpts)
% Code copied from dg_superdeglitch, but instead of removing the glitches
% it just accumulates them on a waveform plot.  Also plots black axis lines
% at +/- <thresh> and at <maxpts>.  Waveforms are shifted vertically so
% that "point a" is always zero.
%INPUTS
% samples: sample data of any size or shape.
% thresh: minimum departure from previous value that qualifies as a glitch.
% maxpts: the maximum number of points that can qualify as a glitch and
%   thus get interpolated away.
%OUTPUT
% numglitch: number of glitches found.
% hF: figure window containing overlaid plots of all glitches found.
%NOTES
%  Based on dg_superdeglitch Rev: 214 (Thu, 26 Mar 2015).

%$Rev: 267 $
%$Date: 2019-09-19 16:17:47 -0400 (Thu, 19 Sep 2019) $
%$Author: dgibson $

samples = reshape(samples,[],1);
difsamp = abs(diff(samples));
hF = figure;
hA = axes('Parent', hF, 'NextPlot', 'add');
idx = 1;
numglitch = 0;
while idx <= (length(samples) - 1)
    % <idx> + 1 points to first point under test as possible glitch.
    if difsamp(idx) <= thresh || isnan(samples(idx))
        % no glitch at <idx>.
        idx = idx + 1;
    else
        % This is a glitch if it's not too long and doesn't end with a
        % NaN. Set <pointa> to point at the putative point "a".
        pointa = idx;
        for k = 1:maxpts
            % <pointb> points to putative point "b":
            pointb = pointa + k + 1;
            if pointb > length(samples) || isnan(samples(pointb))
                break
            end
            foundglitch = false;
            if abs(samples(pointb) - samples(pointa)) <= thresh
                foundglitch = true;
                % found k-point glitch; plot
                numglitch = numglitch + 1;
                endpt = min(pointb + maxpts, length(samples));
                plot( hA, 1 : endpt - pointa + 1, ...
                    samples(pointa : endpt) - samples(pointa) );
                idx = pointb; % so next putative pointa is pointb
                break
            end
        end
        if ~foundglitch
            idx = idx + 1;
        end
    end
end
plot(get(hA, 'XLim'), [1 1] * thresh, 'k');
plot(get(hA, 'XLim'), -[1 1] * thresh, 'k');
plot([1 1] * maxpts, get(hA, 'YLim'), 'k');
title(hA, sprintf('Number of glitches: %d', numglitch));

