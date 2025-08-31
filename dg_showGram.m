function [hI, hCB] = dg_showGram(hDisp, xvals, yvals, data, ...
    titlestr, xlabstr, ylabstr, cbarlabel)
%[hI, hCB] = dg_showGram(hDisp, xvals, yvals, data, ...
%     titlestr, xlabstr, ylabstr, cbarlabel)
% Displays color-scaled <data> with a colorbar in a specified display
% object, which may be either a figure (for single-plot figures) or an axes
% (for subplot figures).  <xvals> and <yvals> are passed on to imagesc,
% which apparently only cares about the first and last and assumes that
% they are linear.  Therefore, a simple check is done on each to see if the
% midpoint is where it should be on a linear scale, and if it isn't, then
% the tick mark labels are created by function 'manticks', which is geared
% primarily towards log scales.  If <cbarlabel> is 'none', then no colorbar
% is created.
%ALSO SEE
% dg_recolorGram

%$Rev: 264 $
%$Date: 2019-07-17 17:41:29 -0400 (Wed, 17 Jul 2019) $
%$Author: dgibson $

if ~isequal([length(yvals), length(xvals)], size(data))
    warning('dg_showGram:axescalib', ...
        'The lengths of the axis values do not match the size of <data>.');
end

switch get(hDisp, 'Type')
    case 'axes'
        hA = hDisp;
    case 'figure'
        hA = axes('Parent', hDisp);
    otherwise
        error('dg_showGram:unknowndisp', ...
            'Unknown display object type: %s', get(hDisp, 'Type'));
end
man_x = ...
    abs(xvals(round(length(xvals)/2)) - (xvals(1) + xvals(end))/2) ...
    > xvals(2)-xvals(1);
man_y = ...
    abs(yvals(round(length(yvals)/2)) - (yvals(1) + yvals(end))/2) ...
    > yvals(2)-yvals(1);
if man_x && man_y
    hI = imagesc(data,'Parent', hA);
    set(hA, 'XTickMode', 'manual');
    set(hA, 'XTickLabelMode', 'manual');
    [t, l] = manticks(xvals);
    set(hA, 'XTick', t);
    set(hA, 'XTickLabel', l);
    set(hA, 'YTickMode', 'manual');
    set(hA, 'YTickLabelMode', 'manual');
    [t, l] = manticks(yvals);
    set(hA, 'YTick', t);
    set(hA, 'YTickLabel', l);
elseif man_x && ~man_y
    hI = imagesc(1:length(xvals), yvals, data,'Parent', hA);
    set(hA, 'XTickMode', 'manual');
    set(hA, 'XTickLabelMode', 'manual');
    [t, l] = manticks(xvals);
    set(hA, 'XTick', t);
    set(hA, 'XTickLabel', l);
elseif ~man_x && man_y
    hI = imagesc(xvals, 1:length(yvals), data,'Parent', hA);
    set(hA, 'YTickMode', 'manual');
    set(hA, 'YTickLabelMode', 'manual');
    [t, l] = manticks(yvals);
    set(hA, 'YTick', t);
    set(hA, 'YTickLabel', l);
else
    hI = imagesc(xvals, yvals, data,'Parent', hA);
end
set(hA, 'YDir', 'normal');
set(get(hA,'XLabel'),'String',xlabstr);
set(get(hA,'XLabel'),'Interpreter','none');
set(get(hA,'YLabel'),'String',ylabstr);
set(get(hA,'YLabel'),'Interpreter','none');
if ~isempty(titlestr)
    title(hA, titlestr, 'Interpreter', 'none');
end
if isequal(cbarlabel, 'none')
    hCB = [];
else
    hCB = colorbar('peer', hA);
    set(get(hCB,'YLabel'),'String',cbarlabel);
    set(get(hCB,'YLabel'),'Interpreter','none');
end
end

function [t, l] = manticks(vals)
% This "should work" to set around ten tick marks at round values of the
% vector <vals>, but it is likely to be pretty only if <vals> is a log
% scale.  The tick coordinates <t> are indices into vals; <l> contains
% label strings.
t = [];
l = {};
if vals(1) <= 0
    % give up - return empties - make no ticks
    return
end
decades = log10(vals(end)/vals(1));
if round(decades) < 1
    t = 1:round(length(vals)/10):length(vals);
    for k = 1:length(t)
        l{k} = sprintf('%g', vals(t(k)));
    end
else
    decPerTick = decades/10;
    scale = 10^floor(log10(vals(1)));
    scaledvals = vals/scale;
    if decPerTick > .5
        decPerTick = round(decPerTick);
        startval = ceil(log10(scaledvals(1)));
        endval = floor(log10(scaledvals(end)));
        tickvals = logspace(startval, endval, ...
            round((endval - startval)/decPerTick) + 1);
    elseif decPerTick > .3
        tickvals = [1 3 10 30 100 300 1000 3000 1e4 3e4 1e5 3e5];
    elseif decPerTick > .17
        tickvals = [1 2 5 10 20 50 100 200 500 1000 2000 5000];
    elseif decPerTick > .1
        tickvals = [1 1.5 2 5 7 10 15 20 50 70 100 150 200 500 700];
    else
        tickvals = [1 1.25 1.5 2 2.5 3 4 5 6 8 10 12.5 15 20 25 30 40 50 60 80];
    end
    % find the smallest tickval that is greater than or equal to vals(1),
    % give or take a factor of a power of ten
    lowtickidx = find(tickvals >= scaledvals(1));
    for tickidx = lowtickidx(1):length(tickvals)
        [m, validx] = min(abs(scaledvals - tickvals(tickidx)));
        if tickvals(tickidx)/scaledvals(end) < 1.001
            t(end+1) = validx;
            l{end+1} = sprintf('%.3g', scale * tickvals(tickidx));
        end
    end
end
end
