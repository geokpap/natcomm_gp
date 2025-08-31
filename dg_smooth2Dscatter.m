function [result, binedgesX, binedgesY] = dg_smooth2Dscatter(coords, values, varargin)
%dg_smooth2Dscatter
% Implements the "Amemori method" of 2D smoothing, to quantize as in a 2D
% histogram, propagate the value into adjacent bins around the point, and
% then compute average values for each bin using only the trial-bins that
% contain values.  Related to 'dg_decisionMatrix' conceptually, but not in
% implementation.  (Not helpful for doing 2D histogramming.)
%INPUTS
% coords: x (col. 1) and y (col. 2) coordinates for a set of observations
%   (one row per observation).
% values: column vector, one value per observation, representing the
%   observation that was associated with the (x, y) coords on the
%   corresponding row of <coords>.
%OUTPUTS
% result: a rectangular numeric array representing a binned, averaged,
%   smoothed rendition of the data in <values> binned according to
%   <coords>.  Bins that contain no <values> are assigned the value <NaN>.
%OPTIONS
% 'bounds', bounds - [minX, maxX, minY, maxY].  Default: [0 1 0 1].
% 'minpts', minpts - the minimum number of data points (rows) to be
%   averaged to produce a non-NaN output value.  Default: 3.
% 'numbins', numbins - number of bins; applies to both axes, so <result> is
%   always square.  Default: 50.
% 'smoothing', smoothing - number of bins by which to propagate each value
%   in all four directions.  Default: 5.
%NOTES
% The calculation is performed via an intermediate 3D array with one square
% plane <numbins> x <numbins> for each observation (row) in <values>.
% Value propagation is done throughout a square of bins centered on the one
% containing the observation.

%$Rev:  $
%$Date:  $
%$Author: dgibson $

minX = 0;
maxX = 1;
minY = 0;
maxY = 1;
minpts = 3;
numbins = 50;
smoothing = 5;

argnum = 0;
while true
    argnum = argnum + 1;
    if argnum > length(varargin)
        break
    end
    if ~ischar(varargin{argnum})
        continue
    end
    switch varargin{argnum}
        case 'bounds'
            argnum = argnum + 1;
            minX = varargin{argnum}(1);
            maxX = varargin{argnum}(2);
            minY = varargin{argnum}(3);
            maxY = varargin{argnum}(4);
        case 'minpts'
            argnum = argnum + 1;
            minpts = varargin{argnum};
        case 'numbins'
            argnum = argnum + 1;
            numbins = varargin{argnum};
        case 'smoothing'
            argnum = argnum + 1;
            smoothing = varargin{argnum};
        otherwise
            error('funcname:badoption', ...
                'The option %s is not recognized.', ...
                dg_thing2str(varargin{argnum}));
    end
end

data = NaN(numbins, numbins, length(values));
binedgesX = linspace(minX, maxX, numbins+1);
binedgesY = linspace(minY, maxY, numbins+1);
for obsidx = 1:length(values)
    Xidx = find(coords(obsidx, 1) >= binedgesX, 1, 'last');
    Yidx = find(coords(obsidx, 2) >= binedgesY, 1, 'last');
    xstart = max(1, Xidx - smoothing);
    xend = min(numbins, Xidx + smoothing);
    ystart = max(1, Yidx - smoothing);
    yend = min(numbins, Yidx + smoothing);
    data(ystart:yend, xstart:xend, obsidx) = values(obsidx);
end
numpts = sum(~isnan(data), 3);
result = mean(data, 3, 'omitnan');
result(numpts < minpts) = NaN;

