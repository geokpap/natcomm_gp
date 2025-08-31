function hF = dg_plotAllClusterFeatures(clusts, features, varargin)
%Plots all unordered pairwise combinations of <features> using
%'plotClusterFeatures'.
%OPTIONS
% 'decimate', N - <N> is the decimation factor for reducing the number of
%   waves in the plot to a manageable number (try to keep it below about
%   1000). If negative, then the decimation factor is calculated separately
%   for each cluster to produce <-N> waveforms.
% 'handle', hF - plot into figure handle <hF> instead of creating figure
%   de novo.
% 'nodisplay' - sets 'Visible' to 'off' when creating figure window.  Has
%   no effect when using pre-existing handle for <hF>.
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
        case 'size'
            argnum = argnum + 1;
            markersize = varargin{argnum};
        otherwise
            error('dg_plotAllClusterFeatures:badoption', ...
                'The option %s is not recognized.', ...
                dg_thing2str(varargin{argnum}));
    end
    argnum = argnum + 1;
end
if isempty(markersize)
    markersize = 6;
end
if isempty(hF)
    hF = figure('Visible', display_str);
end

numcombos = nchoosek(size(features,2), 2);
numrows = ceil(sqrt(numcombos));
plotnum = 0;
if isempty(config)
    sesID = '';
    filename = '';
else
    sesID = config.userparams.sesID;
    filename = config.userparams.basename;
end
for k = 1:size(features,2)-1
    for l = k+1:size(features,2)
        plotnum = plotnum + 1;
        hA = subplot(numrows, numrows, plotnum);
        set(hA, 'Color', 'k');
        dg_plotClusterFeatures(clusts, features(:, [k l]), 'decimate', N, ...
            'handle', hA, 'nolegend', 'size', markersize, 'config', config);
        if isempty(config) || size(features,2) ~= 7
            xlabel(sprintf('feature %d', k));
            ylabel(sprintf('feature %d', l));
        else
            xlabel(config.userparams.featurenames{k});
            ylabel(config.userparams.featurenames{l});
        end
        if plotnum == round(numrows/2)
            title(sprintf('%s %s\n', sesID, filename), ...
                'Interpreter', 'none');
        end
    end
end

