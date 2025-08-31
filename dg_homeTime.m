function [homebase, isathome, xy, splitpoint] = dg_homeTime( ...
    filename, Ts, windowlen, radius, varargin)
%[homebase, isathome, xy, splitpoint] = dg_homeTime( ...
%   filename, Ts, windowlen, radius)
% Reads a whitespace-delimited file containing a file header row, a column
% header header row, and then numeric data in the first seven columns and
% text data in the eighth.  Interprets the column 1 as sample number, col.
% 2 as X, and col. 3 as Y.  Missing sample numbers are interpreted as data
% rows that exactly repeat the previous row.  <Ts> is the sample period in
% seconds.  <windowlen> is the number of samples to use for calculating
% individual "home base time" values, which are expressed as the fraction
% of samples where the position is at home base. <radius> is the radius
% within which the position is considered to be "at home base".  The
% location of "home base" is defined as the most frequent (X, Y) pair in
% the file.  If there are several such pairs that are within <radius> of
% each other, then "home base" is the average position of all such pairs.
% Scores each sample as being within <radius> of home base or not and plots
% a <windowlen>-point trailing moving average.
% OUTPUTS
%   homebase - a 1x2 array containing x,y coordinates of home base
%   isathome - a column vector of length equal to number of data points
%       (after expanding repeated points) that is true for points where the
%       animal was within <radius> of home base.
%   xy - an Nx2 array of x,y coordinates after expanding repeated points.
%   splitpoint - the row of xy that corresponds to splitT if 'splitT' is
%       given; otherwise, [].
% OPTIONS:
%   'axes', hA - hA is a two-column array of axes handles.  Plots
%       2D occupancy histogram in hA(:,1), dg_plotHomeTime in hA(:,2).  hA
%       should have one row if 'spliT' is not specified, two rows if
%       'spliT' is specified.
%   'label', labelval - if <labelval> is 1, puts title and full labels on
%       plots. If 0, puts no labels at all.  If a string, labels the y-axis
%       only with the string.
%   'splitT', splitT - <splitT> is a time specified in minutes that is used
%       to divide the data into a "before" period and an "after" period.
%       The entire analysis is repeated for both periods, and the output
%       values <homebase> and <isathome> are replaced by cell vectors,
%       where each cell contains the value for the corresponding period.
%   'xmax', xmax - sets x-axis scale.
%   'ymax', ymax - sets y-axis scale.
% NOTES:
% (1) The file header row (which seems to be the filename without the
% extension) must NOT contain any whitespace characters.
% (2) For comparisons between <radius> and distances that are calculated as
% floating-point, (radius + 1e-6) is used instead of radius to guarantee
% that even if there are truncation errors, distances that are exactly
% equal to radius will still be counted as part of the home base.

%$Rev: 24 $
%$Date: 2009-03-31 21:51:08 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

hA = [];
imageaxes = [];
labelval = 1;
lineplotaxes = [];
splitT = [];
xmax= [];
ymax = [];
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'axes'
            argnum = argnum + 1;
            hA = varargin{argnum};
            if length(hA) ~= 2 || ~isa(hA, 'double')
                error('dg_homeTime:hA', ...
                    'hA must be a two-element vector of axes handles.' );
            end
            imageaxes = hA(:,1);
            lineplotaxes = hA(:,2);
        case 'label'
            argnum = argnum + 1;
            labelval = varargin{argnum};
        case 'splitT'
            argnum = argnum + 1;
            splitT = varargin{argnum};
            if ~isequal(size(splitT), [1 1]) || ~isnumeric(splitT)
                error('dg_homeTime:splitT', ...
                    'splitT value should be a numeric scalar: %s', ...
                    dg_thing2str(splitT) );
            end
        case 'xmax'
            argnum = argnum + 1;
            xmax = varargin{argnum};
            if ~isequal(size(xmax), [1 1]) || ~isnumeric(xmax)
                error('dg_homeTime:xmax', ...
                    'xmax value should be a numeric scalar: %s', ...
                    dg_thing2str(xmax) );
            end
        case 'ymax'
            argnum = argnum + 1;
            ymax = varargin{argnum};
            if ~isequal(size(ymax), [1 1]) || ~isnumeric(ymax)
                error('dg_homeTime:ymax', ...
                    'ymax value should be a numeric scalar: %s', ...
                    dg_thing2str(ymax) );
            end
        otherwise
            error('dg_homeTime:badoption', ...
                ['The option "' dg_thing2str(varargin{argnum}) '" is not recognized.'] );
    end
    argnum = argnum + 1;
end
splitpoint = round((splitT * 60) / Ts);

fid = fopen(filename);
if fid == -1
    error('dg_homeTime:file', ...
        'Could not open %s', filename);
end
C = textscan(fid, '%s');

% Unpack the relevant portion of the file and convert to numeric array
C = C{1};
C(1:9) = [];
C = reshape(C, 8, [])';
nums = cell2mat(dg_mapfunc(@str2num, C(:, 1:3)));

% Find sequence gaps, and create array of (x,y) coords with gaps filled in.
% I define the first element not to be a gap.
gapidx = find([0; nums(2:end,1) ~= nums(1:end-1,1) + 1]);
xy(nums(:,1), :) = nums(:, 2:3);
for gaprow = reshape(gapidx, 1, [])
    gap = (nums(gaprow-1, 1) + 1) : (nums(gaprow, 1) - 1);
    xy(gap, :) ...
        = repmat(nums(gaprow-1, 2:3), length(gap), 1);
end

[pathstr, name] = fileparts(filename);
if isempty(splitT)
    homebase = findHomeBase(xy, name, radius, imageaxes, ...
        labelval, xmax, ymax);
    isathome = dg_plotHomeTime(homebase, xy, name, Ts, windowlen, ...
        radius, 1, lineplotaxes, labelval);
else
    homebase{1} = findHomeBase(xy(1:splitpoint, :), ...
        sprintf('%s %.0f-%.0f s', name, 0, (splitpoint-1)*Ts), ...
        radius, imageaxes(1), labelval, xmax, ymax);
    homebase{2} = findHomeBase(xy(splitpoint+1:end, :), ...
        sprintf('%s %.0f-%.0f s', name, splitpoint*Ts, (size(xy,1)-1)*Ts), ...
        radius, imageaxes(2), labelval, xmax, ymax);
    isathome{1} = dg_plotHomeTime(homebase{1}, xy(1:splitpoint, :), ...
        sprintf('%s %.0f-%.0f s', name, 0, (splitpoint-1)*Ts), ...
        Ts, windowlen, radius, 1, lineplotaxes(1), labelval);
    isathome{2} = dg_plotHomeTime(homebase{2}, xy(splitpoint+1:end, :), ...
        sprintf('%s %.0f-%.0f s', name, splitpoint*Ts, (size(xy,1)-1)*Ts), ...
        Ts, windowlen, radius, splitpoint+1, lineplotaxes(2), labelval);
end

end


function homebase = findHomeBase(xy, name, radius, imageaxes, ...
    labelval, xmax, ymax)
% Make 2-D occupancy histogram by converting coordinates into array
% indices and then accumulating counts. Zero might be a legitimate
% coordinate.  
if isempty(xmax)
    xmax = max(xy(:,1));
elseif any(xy(:,1)>xmax)
    error('dg_homeTime:xmax', ...
        'There are x values greater than xmax');
end
if isempty(ymax)
    ymax = max(xy(:,2));
elseif any(xy(:,2)>ymax)
    error('dg_homeTime:ymax', ...
        'There are y values greater than ymax');
end
xidx = 2 * xy(:,1) + 1;
yidx = 2 * xy(:,2) + 1;
xyhist = zeros(2 * ymax + 1, 2 * xmax + 1);
for point = 1:size(xy, 1)
    xyhist(yidx(point), xidx(point)) = ...
        xyhist(yidx(point), xidx(point)) + 1;
end

% Plot occupancy
if isempty(imageaxes)
    figure;
    imageaxes = axes;
    titlestr = sprintf('Space Occupancy Histogram %s', name);
else
    titlestr = name;
end
imagesc(0:xmax, 0:ymax, log10(xyhist), ...
    'Parent', imageaxes);
cbaraxes = colorbar('peer', imageaxes);
switch labelval
    case 1
        title(imageaxes, titlestr);
        ylabel(cbaraxes, 'Log Number of Tracker Samples');
    case 0
        set(imageaxes, 'XTickLabel', []);
    otherwise
        set(imageaxes, 'XTickLabel', []);
        ylabel(imageaxes, labelval);
end

% Find home base
maxcount = max(xyhist(:));
[homebase(:,2), homebase(:,1)] = find(xyhist==maxcount);
if size(homebase, 1) > 1
    avghome = mean(homebase);
    if all(sqrt( (homebase(:,1) - avghome(1)).^2 ...
            + (homebase(:,2) - avghome(2)).^2 ) <= (radius + 1e-6))
        homebase = avghome;
    else
        error('dg_homeTime:homes', ...
            'There are multiple home bases.' );
    end
end
homebase = (homebase - [1 1])/2;

end
