function hL = lfp_plotTMazeEvents(hA, evtIDs, ...
    xposn, yposn, outliersx, outliersy, legendflag)
% Plots color-coded event markers on tracker plot, using different symbols
% to distinguish outliers. 
%INPUTS
% hA - axes handle into which to plot.
% evtIDs - a cell vector of lists of functionally equivalent event IDs.
% xposn, yposn, outliersx, outliersy - coordinates of event markers,
%   separated into outliers and non-outliers.  These are 2-D
%   arrays in (evtID x trial) format.  All four arrays are the same size
%   and contain NaNs where are no values to plot.
% legendflag - set true to include a legend.  Gets very confused if there
%   is anything more than tracker data in the plot, or even if the tracker
%   data were plotted as more than one single line object.
%OUTPUT
% hL - a column vector of line handles for all the line objects plotted by
%   this function (one line object for each element of <evtIDs>).

%$Rev: 85 $
%$Date: 2009-10-10 00:57:03 -0400 (Sat, 10 Oct 2009) $
%$Author: dgibson $

lfp_declareGlobals;

if length(evtIDs) ~= size(xposn,1)
    error('lfp_plotTMazeEvents:evtIDs', ...
        'Position arrays must have one row per evtID');
end
if ~isequal(size(xposn), size(yposn)) || ...
        ~isequal(size(xposn), size(outliersx)) || ...
        ~isequal(size(xposn), size(outliersy))
    error('lfp_plotTMazeEvents:posns', ...
        'All four position arrays must be the same size');
end

set(hA, 'NextPlot', 'add');
hL = [];

c = get(hA, 'ColorOrder');
blackrow = find(c(:,1)<.3 & c(:,2)<.3 & c(:,3)<.3);
if ~isempty(blackrow)
    c(blackrow,:) = [0.7 0.7 0.7];
end;
legendstr = {'tracker'};
for evtIDidx = 1:length(evtIDs)
    if ~all(isnan(xposn(evtIDidx, :)))
        hL(end+1,1) = plot(hA, xposn(evtIDidx, :), yposn(evtIDidx, :), ...
            'Color', c(mod(evtIDidx-1, size(c,1)) + 1, :), ...
            'Marker', '.', 'LineStyle', 'none', ...
            'MarkerSize', 18);
        legendstr{end+1} = sprintf('evt %s', ...
            dg_thing2str(evtIDs{evtIDidx}));
    end
end
if legendflag
    legend(legendstr, 'Location', 'NorthWest');
end
% Plot outliers separately but with same colors so they don't need
% separate legend entries:
for evtIDidx = 1:length(evtIDs)
    if ~all(isnan(outliersx(evtIDidx, :)))
        hL(end+1,1) = plot(hA, outliersx(evtIDidx, :), ...
            outliersy(evtIDidx, :), ...
            'Color', c(mod(evtIDidx-1, size(c,1)) + 1, :), ...
            'Marker', '*', 'LineStyle', 'none' );
    end
end