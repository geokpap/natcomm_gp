function hF = lfp_plot_hht(S_mean, timebinctrs, freqs, plotdata)
%hF = lfp_plot_hht(S_mean, timebinctrs, freqs, plotdata)
% All inputs are as returned by lfp_hht.

%$Rev: 147 $
%$Date: 2010-07-27 15:23:57 -0400 (Tue, 27 Jul 2010) $
%$Author: dgibson $

hF = figure;
if plotdata.dbflag
    cscalestr = 'Power, dB (pre)';
elseif plotdata.dbpostflag
    cscalestr = 'Power, dB (post)';
elseif plotdata.dbprepostflag
    cscalestr = 'Power, dB (prepost)';
else
    cscalestr = 'Power (linear)';
end
[hI, hCB] = dg_showGram(hF, timebinctrs, freqs, S_mean, ...
    '', 'Time, s', 'Frequency, Hz', cscalestr);
if isempty(plotdata.autocolor)
    if ~isempty(plotdata.clim)
        dg_recolorGram(hCB, plotdata.clim, hI);
    end
else
    clim = [-plotdata.autocolor 0] + max(S_mean(:));
    dg_recolorGram(hCB, clim, hI);
end
hA = get(hI, 'Parent');
lfp_createFigTitle(hA, [],[],[],[],[], 'figtitle', plotdata.figtitle);

if ~isempty(plotdata.semlevel)
    caxis(hA, caxis); % pre-de-botch the color scale
    xi = linspace(timebinctrs(1), timebinctrs(end), length(timebinctrs)*4);
    yi = linspace(freqs(1), freqs(end), length(freqs)*4)';
    semi = interp2(reshape(timebinctrs, 1, []), ...
        reshape(freqs, [], 1), plotdata.S_std/sqrt(plotdata.N), xi, yi);
    set(hA, 'NextPlot', 'add');
    if numel(plotdata.semlevel) > 1
        graylevel = linspace(1, .3, numel(plotdata.semlevel));
    else
        graylevel = 1;
    end
    for levnum = 1:numel(plotdata.semlevel)
        contour(hA, xi, yi, semi, [1 1] * plotdata.semlevel(levnum), ...
            'Color', [1 1 1] * graylevel(levnum));
    end
end
if length(plotdata.trials) == 1
    lfp_plotEvtMarkers( hA, ...
        plotdata.trialevts(plotdata.evts2mark,:), ...
        'reftime', plotdata.reftime, ...
        'bigevts', plotdata.bigevts );
else
    lfp_plotEvtMarkers(hA, [0 lfp_AlignmentRef(1)]);
end

