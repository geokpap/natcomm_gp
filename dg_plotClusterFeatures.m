function hF = dg_plotClusterFeatures(clusts, features, varargin)
%INPUTS
% features: one row per observation for some set of observations, one
%   feature per column, only first two columns are used.
% clusts: cell vector of arrays of numeric indices into <feature>.
%NOTES
% Plots points that do not belong to any cluster as white.
%OPTIONS
% 'decimate', N - <N> is the decimation factor for reducing the number of
%   waves in the plot to a manageable number (try to keep it below about
%   1000). If negative, then the decimation factor is calculated separately
%   for each cluster to produce <-N> waveforms.
% 'handle', hF - plot into figure or axes handle <hF> instead of creating
%   figure and axes de novo.
% 'nodisplay' - sets 'Visible' to 'off' when creating figure window.  Has
%   no effect when using pre-existing handle for <hF>.
% 'nolegend' - suppresses display of plot legend and title.
% 'size', markersize - points are normally plotted at Matlab default size.
%   This overrides that.  <markersize> can be a number, or 'various'; in
%   the latter case, markers start out large and get progressively smaller
%   until the last set of points is plotted.

%$Rev: 301 $
%$Date: 2022-09-22 18:56:54 -0400 (Thu, 22 Sep 2022) $
%$Author: dgibson $

config = [];
markersize = [];
display_str = 'on';
hF = [];
N = 1;
plotlegend = true;
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'config'
            argnum = argnum + 1;
            config = varargin{argnum};
        case 'decimate'
            argnum = argnum + 1;
            N = varargin{argnum};
        case 'handle'
            argnum = argnum + 1;
            hF = varargin{argnum};
        case 'nodisplay'
            display_str = 'off';
        case 'nolegend'
            plotlegend = false;
        case 'size'
            argnum = argnum + 1;
            markersize = varargin{argnum};
        otherwise
            error('dg_plotClusterFeatures:badoption', ...
                'The option %s is not recognized.', ...
                dg_thing2str(varargin{argnum}));
    end
    argnum = argnum + 1;
end
defaultsize = 6;
if isequal(markersize, 'various')
    markersize = linspace( defaultsize + 6*length(clusts), ...
        defaultsize, length(clusts)+1 );
elseif isempty(markersize)
    markersize = defaultsize * ones(length(clusts)+1, 1);
else
    markersize = markersize * ones(length(clusts)+1, 1);
end

if isempty(hF)
    hF = figure('Visible', display_str);
    hA = axes('Parent', hF, 'NextPlot', 'add', 'Color', 'k');
else
    if isequal(get(hF, 'Type'), 'axes')
        hA = hF;
        set(hA, 'NextPlot', 'add');
    else
        hA = axes('Parent', hF, 'NextPlot', 'add', 'Color', 'k');
    end
end
if ~iscell(clusts)
    error('dg_plotClusterFeatures:clusts', ...
        '<clusts> must be a cell vector of arrays of numeric indices.');
end

allclusts = [];
for k = 1:length(clusts)
    if ~isnumeric(clusts{k})
        error('dg_plotClusterFeatures:notnum', ...
            'Cluster %d is not numeric.', k);
    end
    allclusts = [ allclusts
        clusts{k}(:)
        ]; %#ok<*AGROW>
end
clust0 = setdiff(1:size(features,1), allclusts);
colors = get(hA, 'ColorOrder');
% Set Plexon OFS colors:
colors(1:5, :) = [
    255 255 0
    0 255 0
    0 255 255
    255 0 0
    192 192 0
    ] / 255;
% color0 = [192 192 192] / 255; % exact match to Plexon
% under "C" in Matlab color palette:
color0 = [0.831372549019608 0.815686274509804 0.784313725490196];
hLclust = [];
if isempty(config)
    sesID = '';
    filename = '';
else
    sesID = config.userparams.sesID;
    filename = config.userparams.basename;
end
for clustidx = 0:length(clusts)
    if clustidx == 0
        obsidx = clust0;
        rgb = color0;
    else
        obsidx = clusts{clustidx};
        coloridx = mod(clustidx-1, size(colors,1))+1;
        rgb = colors(coloridx,:);
    end
    if N < 0
        incr = round(-numel(obsidx)/N);
    else
        incr = N;
    end
    if isempty(markersize)
        opts = {};
    else
        opts = {'MarkerSize', markersize(clustidx+1)};
    end
    hL = plot(hA, ...
        features(obsidx(1:incr:end),1), features(obsidx(1:incr:end),2), ...
        '.', 'color', rgb, opts{:});
    if clustidx > 0
        hLclust(clustidx) = hL(1);
        clustIDvalstr{clustidx} = int2str(clustidx);
    end
end
if plotlegend
    title(sprintf('%s %s\n', sesID, filename), ...
                'Interpreter', 'none');
    if ~isempty(hLclust)
        legend(hLclust, clustIDvalstr, 'TextColor', 'w');
    end
end


