function hF = dg_plotWavePrctiles(timepts, wavedata, varargin)
% Plots prctiles [0 5 50 95 100] across all spikes.
%INPUTS
% timepts: vector of time points corresponding to the points in
%   <wavedata>.
% wavedata: rectangular array of numeric data with one row for each
%	observation, one column per sample.
%OUTPUTS
% hF: handle to figure plotted.
%OPTIONS
% 'color', rgb - plots all lines using color specified by RGB triple <rgb>.
% 'handle', hF - plot into figure or axes handle <hF> instead of creating
%   figure and axes de novo.
% 'nodisplay' - sets 'Visible' to 'off' when creating figure window.  Has
%   no effect when using pre-existing handle for <hF>.
%NOTES
% Derived from 'dg_plotClusterWaves'.

%$Rev: 309 $
%$Date: 2024-02-01 14:43:52 -0500 (Thu, 01 Feb 2024) $
%$Author: dgibson $

display_str = 'on';
hF = [];
argnum = 1;
prctiles = [0 5 50 95 100];
rgb = [0 0.447 0.741];
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'color'
            argnum = argnum + 1;
            rgb = varargin{argnum};
        case 'handle'
            argnum = argnum + 1;
            hF = varargin{argnum};
        case 'nodisplay'
            display_str = 'off';
        otherwise
            error('dg_plotWavePrctiles:badoption', ...
                'The option %s is not recognized.', ...
                dg_thing2str(varargin{argnum}));
    end
    argnum = argnum + 1;
end

if ~isequal(size(rgb), [1 3])
    error('dg_plotWavePrctiles:rgb', ...
        '<rgb> must be a 3 element row vector.');
end

if isempty(hF)
    hF = figure('Visible', display_str);
    hA = axes('Parent', hF, 'NextPlot', 'add');
else
    if isequal(get(hF, 'Type'), 'axes')
        hA = hF;
        set(hA, 'NextPlot', 'add');
    else
        hA = axes('Parent', hF, 'NextPlot', 'add', 'Color', 'k');
    end
end
waves = prctile(wavedata, prctiles);
medidx = 3;
plot(hA, timepts, waves(medidx,:)', ...
    'color', rgb, ...
    'LineWidth', 2); 
plot(hA, timepts, waves([2 4],:)', 'color', rgb, 'LineWidth', 1);
plot(hA, timepts, waves([1 5],:)', 'color', rgb, 'LineWidth', 0.5);

