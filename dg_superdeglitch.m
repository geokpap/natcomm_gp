function samples = dg_superdeglitch(samples, thresh, maxpts)
%samples = dg_superdeglitch(samples, thresh, maxpts)
%INPUTS
% samples: sample data of any size or shape.
% thresh: minimum departure from previous value that qualifies as a glitch.
% maxpts: the maximum number of points that can qualify as a glitch and
%   thus get interpolated away.
%OUTPUT
%NOTES
%  Returns the same sample data that was passed in, modified by linearly
%  interpolating glitches.  A k-point glitch is defined as a series of k
%  points that differ from a preceding point "a" by more than <thresh>, but
%  which is followed by a point "b" that differs from the predecessor point
%  "a" by less than or equal to <thresh>. It is explicitly NOT required
%  that the k members of the k-point glitch be within <thresh> of each
%  other, only that they be more than <thresh> from "a".  This makes it
%  possible to remove random noise bursts.  The interpolation extends from
%  the last point before the glitch through (including) the first point
%  after the glitch, so all glitch points get replaced by the
%  interpolation.
%   NaN values do not qualify as glitches, but they also don't
%  qualify as endpoints of a sequence to interpolate, and will prevent any
%  glitches that contain them from being removed.
%   This version has been completely rewritten to use less memory.
%  However, it does still have some memory issues, to wit:  the amount of
%  additional memory it requires is about twice the size of <samples>,
%  which means that the amount of memory allocated to the Matlab running it
%  will be at least three times the size of <samples> plus the size of a
%  freshly started Matlab (which is usually around 2 GB).

%$Rev: 267 $ $Date: 2019-09-19 16:17:47 -0400 (Thu, 19 Sep 2019) $ $Author:
%dgibson $

origsize = size(samples);
samples = reshape(samples,[],1);
difsamp = abs(diff(samples));
idx = 1;
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
                % found k-point glitch; interpolate
                samples(pointa:pointb) = linspace(samples(pointa), ...
                    samples(pointb), pointb - pointa + 1);
                idx = pointb; % so next putative pointa is pointb
                break
            end
        end
        if ~foundglitch
            idx = idx + 1;
        end
    end
end

samples = reshape(samples, origsize);

