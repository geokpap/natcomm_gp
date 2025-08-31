function hF = dg_plotClusterWavesOverlay(clustIDs, wavedata, varargin)
%INPUTS
% clustIDs: may be numeric, in which case it is spike-by-spike cluster IDs;
%   or may be a cell vector, in which case it is a column vector each of
%   whose cells contains a vector of indices into <wavedata>.  When using
%   the cell vector format, it is possible for a waveform to belong to more
%   than cluster, in which case higher cluster numbers overwrite lower
%   cluster numbers, and a warning is raised.  Any waveforms that do not
%   belong to any cluster in the cell vector format are designated as
%   cluster "0".
% wavedata: rectangular array of numeric data with one row for each
%	observation, one column per sample.
%OPTIONS
% 'decimate', N - <N> is the decimation factor for reducing the number of
%   waves in the plot to a manageable number (try to keep it below about
%   1000). If negative, then the decimation factor is calculated separately
%   for each cluster to produce <-N> waveforms.
% 'handle', hF - plot into figure or axes handle <hF> instead of creating
%   figure and axes de novo.
% 'nodisplay' - sets 'Visible' to 'off' when creating figure window.  Has
%   no effect when using pre-existing handle for <hF>.
%NOTES
% Plots waveforms for cluster 0 as gray.

%$Rev: 301 $
%$Date: 2022-09-22 18:56:54 -0400 (Thu, 22 Sep 2022) $
%$Author: dgibson $

config = [];
display_str = 'on';
hF = [];
N = 1;
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
        otherwise
            error('dg_plotClusterWavesOverlay:badoption', ...
                'The option %s is not recognized.', ...
                dg_thing2str(varargin{argnum}));
    end
    argnum = argnum + 1;
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

if iscell(clustIDs)
    % convert <clustIDs> format to numeric:
    newclustIDs = zeros(size(wavedata,1),1);
    for idx = 1:length(clustIDs)
        if any(newclustIDs(clustIDs{idx}))
            warning('dg_plotClusterWavesOverlay:clustIDs', ...
                '<clustIDs> contains overlapping clusters.');
        end
        newclustIDs(clustIDs{idx}) = idx;
    end
    clustIDs = newclustIDs;
elseif ~isnumeric(clustIDs)
    error('dg_plotClusterWavesOverlay:clustIDs', ...
        'bad data type for <clustIDs>');
end
clustIDvals = setdiff(unique(clustIDs), 0);
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
titlestr = 'Spike Counts';
gotclust0 = 1;
numpts = size(wavedata, 2);
if isempty(config)
    xvals = 1:numpts;
    sesID = '';
    filename = '';
else
    sesID = config.userparams.sesID;
    filename = config.userparams.basename;
    % time values in microseconds
    wavdur = 1e6 * (config.userparams.wavelen - 1) / ...
        config.userparams.samplefreq;
    sampleperiod = wavdur / (numpts - 1);
    ovrsampalign = config.userparams.align * numpts / ...
        config.userparams.wavelen;
    xvals = (-(ovrsampalign-1) : (numpts - ovrsampalign)) * sampleperiod;
    % voltage values in microvolts
    wavedata = 1e6 * config.userparams.bitvolts * wavedata;
end
for clustidx = 0:length(clustIDvals)
    if clustidx == 0
        rgb = color0;
        clustid = 0;
    else
        coloridx = mod(clustidx-1, size(colors,1))+1;
        rgb = colors(coloridx,:);
        clustid = clustIDvals(clustidx);
    end
        clustwaveidx = find(clustIDs == clustid);
    titlestr = sprintf( '%s / %d:%d', titlestr, ...
        clustid, length(clustwaveidx) );
    if N < 0
        incr = max(1, round(-length(clustwaveidx)/N));
    else
        incr = N;
    end
    waves = wavedata(clustwaveidx(1:incr:end), :);
    if clustid == 0 && isempty(waves)
        gotclust0 = 0;
    else
        hL = plot(hA, xvals, waves', 'color',rgb);
        hLclust(clustidx + gotclust0) = hL(1); %#ok<AGROW>
        clustIDvalstr{clustidx + gotclust0} = int2str(clustid); %#ok<AGROW>
    end
end
legend(hA, hLclust, clustIDvalstr, 'TextColor', 'w');
title(sprintf('%s %s\n%s', sesID, filename, titlestr), ...
    'Interpreter', 'none');
if isempty(config)
    ylabel('AD units');
    xlabel('samples');
else
    % set X scale
    xlabel('microseconds');
    set(hA, 'XLim', xvals([1 end]));
    % set Y scale
    ylabel('microvolts');
    set(hA, 'YLim', [-1 1] * 1e6 * config.userparams.maxval ...
        * config.userparams.bitvolts );
    % draw trigger threshold
    plot(get(hA, 'XLim'), [1 1] * config.userparams.thresh, 'Color', 'w');
    % draw time zero
    plot([0 0], get(hA, 'YLim'), 'Color', 'w');
end


