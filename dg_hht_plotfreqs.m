function hF = dg_hht_plotfreqs(f, A, Ts, startpt, figtitle, flim, maxA, ...
    flags, plotimf)
% Create figure containing plot of instantaneous frequencies vs. time for
% dg_hht.
%INPUTS
% f: frequency in radians per sample in IMFs x samples format.
% A: amplitude
% Ts: sample period
% startpt: position of first sample point in plot relative to nominal time
%   0, in samples.
% figtitle: indeed
% flim: frequency limits for plots; points outside this range are not
%   plotted.  Optional; default is the range of values in <f>.
% maxA: value of <A> that corresponds to maximum marker size in the
%   freq-amplitude traces.  Optional; default is to make the highest value
%   of <A> in the plot be maximum size.
% flags: first element enables line plot, second element enables
%   freq-amplitude plot.  Optional; default is to enable both.
% plotimf: logical vector to individually enable plotting of selected IMFs.
%   Optional; default is to plot all.
%OUTPUTS
% hF: vector of handles to figs created, hF(1) is freqs, hF(2) is freq-amp.

%$Rev: 65 $
%$Date: 2010-08-03 20:36:46 -0400 (Tue, 03 Aug 2010) $
%$Author: dgibson $

f = f / (Ts*2*pi); % f now in Hz

if nargin < 9
    plotimf = true(size(f,1), 1);
end

if nargin < 8
    flags = [true true];
end

if nargin < 7
    maxA = [];
end

if nargin < 6
    flim = [];
end

timepts = Ts * ((1:size(f,2)) - 1 + startpt);
if isempty(flim)
    flim = [min(f(:)) max(f(:))];
end
finflim = f>=flim(1) | f<=flim(2);

if flags(1)
    % freq traces:
    hF(1) = figure;
    hA = axes('Parent', hF(1), 'NextPlot', 'add');
    colors = get(hA,'ColorOrder');
    f(~finflim) = NaN;
    for cidx = 1:size(f,1)
        if plotimf(cidx)
            mycolor = colors(mod(cidx-1, size(colors,1)) + 1, :);
            hL = plot(hA, timepts, f(cidx,:)', 'Color', mycolor);
            set(hL, 'ButtonDownFcn', sprintf('disp(''IMF#%d'');', cidx));
        end
    end
    xlim(timepts([1 end]));
    ylim(flim);
    grid(hA, 'on');
    title(hA, figtitle, 'Interpreter', 'none');
    xlabel('Time, s');
    ylabel('Frequency, Hz');
    legend(hA, cellstr(num2str(find(plotimf))));
end
if flags(2)
    % freq-amplitude traces:
    hF(2) = figure;
    hA = axes('Parent', hF(2), 'NextPlot', 'add');
    colors = get(hA,'ColorOrder');
    if isempty(maxA)
        maxA = max(A(finflim));
    end
    for cidx = 1:size(f,1)
        if plotimf(cidx)
            markersz = round(100 * A(cidx,:) / maxA);
            markersz(markersz<=0) = 1;
            scatter( timepts, ...
                f(cidx,:), markersz, ...
                colors(mod(cidx-1, size(colors,1)) + 1, :), 'filled' );
        end
    end
    xlim(timepts([1 end]));
    ylim(flim);
    grid(hA, 'on');
    title(hA, sprintf('%s maxA=%.2g', figtitle, maxA), 'Interpreter', 'none');
    xlabel('Time, s');
    ylabel('Frequency, Hz');
    legend(hA, cellstr(num2str(find(plotimf))));
end
