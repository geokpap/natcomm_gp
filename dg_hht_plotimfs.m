function hF = dg_hht_plotimfs(imf, Ts, startpt, figtitle, separate, ...
    plotimf)
% Create figure containing plot of IMF waveforms vs. time for dg_hht.
%OUTPUT
% hF: a vector of handles to the figures created

%$Rev: 65 $
%$Date: 2010-08-03 20:36:46 -0400 (Tue, 03 Aug 2010) $
%$Author: dgibson $

timepts = Ts * ((1:size(imf,2)) - 1 + startpt);

if nargin < 6
    plotimf = true(size(imf,1), 1);
end

numimfs = sum(plotimf);
imfs2plot = find(plotimf);

if separate
    plotsperfig = 4;
    numfigs = (ceil(numimfs/plotsperfig));
    for fignum = 1:numfigs
        if fignum < numfigs
            numsubplots = plotsperfig;
        else
            numsubplots = mod(numimfs, plotsperfig);
            if numsubplots == 0
                numsubplots = plotsperfig;
            end
        end
        hF(fignum) = figure; %#ok<AGROW>
        for plotnum = 1:numsubplots
            imfnum = imfs2plot((fignum-1) * plotsperfig + plotnum);
            hA = subplot(plotsperfig, 1, plotnum, 'Parent', hF(fignum));
            hL = plot(hA, timepts, imf(imfnum,:));
            set(hL, 'ButtonDownFcn', sprintf('disp(''IMF#%d'');', imfnum));
            xlim(timepts([1 end]));
            grid(hA, 'on');
            if plotnum == numsubplots
                xlabel('Time, s');
            end
            ylabel(sprintf('IMF#%d', imfnum));
        end
        hA = subplot(plotsperfig, 1, 1);
        title(hA, figtitle, 'Interpreter', 'none');
    end
else
    hF = figure;
    hA = axes('Parent', hF(1), 'NextPlot', 'add');
    colors = get(hA,'ColorOrder');
    for cidx = 1:size(imf,1)
        if plotimf(cidx)
            mycolor = colors(mod(cidx-1, size(colors,1)) + 1, :);
            hL = plot(hA, timepts, imf(cidx,:)', 'Color', mycolor);
            set(hL, 'ButtonDownFcn', sprintf('disp(''IMF#%d'');', cidx));
        end
    end
    xlim(timepts([1 end]));
    grid(hA, 'on');
    title(hA, figtitle, 'Interpreter', 'none');
    xlabel('Time, s');
    ylabel('IMFs');
    legend(hA, cellstr(num2str(imfs2plot)));
end
