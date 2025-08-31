function hF = dg_plotClusterWaves(clustIDs, wavedata, varargin)
% Plots prctiles [0 5 50 95 100] across all spikes.
%INPUTS
% clustIDs: may be numeric, in which case it is spike-by-spike cluster IDs;
%   or may be a cell vector, in which case it is a column vector each of
%   whose cells contains a vector of indices into <wavedata>.  When using
%   the cell vector format, it is possible for a spike to belong to more
%   than cluster, in which case higher cluster numbers overwrite lower
%   cluster numbers, and a warning is raised.  Any waveforms that do not
%   belong to any cluster in the cell vector format are designated as
%   cluster "0".
% wavedata: rectangular array of numeric data with one row for each
%	observation, one column per sample.
%OUTPUTS
% hF: handle to figure plotted.
%OPTIONS
% 'filename', filename - <filename> is shown in the title.
% 'handle', hF - plot into figure or axes handle <hF> instead of creating
%   figure and axes de novo.
% 'nodisplay' - sets 'Visible' to 'off' when creating figure window.  Has
%   no effect when using pre-existing handle for <hF>.
% 'session', sesID - <sesID> is shown in the title.
%NOTES
% See 'dg_plotWavePrctiles' for a more general-purpose prctile plotter.

%$Rev: 308 $
%$Date: 2024-02-01 14:40:45 -0500 (Thu, 01 Feb 2024) $
%$Author: dgibson $

config = [];
display_str = 'on';
filename = '';
hF = [];
sesID = '';
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'config'
            argnum = argnum + 1;
            config = varargin{argnum};
        case 'filename'
            argnum = argnum + 1;
            filename = varargin{argnum};
        case 'handle'
            argnum = argnum + 1;
            hF = varargin{argnum};
        case 'nodisplay'
            display_str = 'off';
        case 'session'
            argnum = argnum + 1;
            sesID = varargin{argnum};
        otherwise
            error('dg_plotClusterWaves:badoption', ...
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
        set(hA, 'NextPlot', 'add', 'Color', 'k');
    else
        hA = axes('Parent', hF, 'NextPlot', 'add', 'Color', 'k');
    end
end
if iscell(clustIDs)
    % convert <clustIDs> format to numeric:
    newclustIDs = zeros(size(wavedata,1),1);
    for idx = 1:length(clustIDs)
        if any(newclustIDs(clustIDs{idx}))
            warning('dg_plotClusterWaves:clustIDs', ...
                '<clustIDs> contains overlapping clusters.');
        end
        newclustIDs(clustIDs{idx}) = idx;
    end
    clustIDs = newclustIDs;
elseif ~isnumeric(clustIDs)
    error('dg_plotClusterWaves:clustIDs', ...
        'bad dat type for <clustIDs>');
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
numpts = size(wavedata, 2);
if isempty(config)
    xvals = 1:numpts;
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
        clustid = 0;
        isinclust = clustIDs == 0;
        rgb = color0;
    else
        coloridx = mod(clustidx-1, size(colors,1))+1;
        clustid = clustIDvals(clustidx);
        isinclust = clustIDs == clustIDvals(clustidx);
        rgb = colors(coloridx,:);
    end
    waves = prctile(wavedata(isinclust,:), [0 5 50 95 100]);
    hLmed(clustidx+1) = plot(hA, xvals, waves(3,:)', ...
        'color', rgb, ...
        'LineWidth', 2); %#ok<AGROW>
    plot(hA, xvals, waves([2 4],:)', 'color', rgb, 'LineWidth', 1);
    plot(hA, xvals, waves([1 5],:)', 'color', rgb, 'LineWidth', 0.5);
    clustIDvalstr{clustidx+1} = int2str(clustid); %#ok<AGROW>
end
title(sprintf('%s %s', sesID, filename), 'Interpreter', 'none');
legend(hLmed, clustIDvalstr, 'TextColor', 'w', 'Color', 'k');
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

