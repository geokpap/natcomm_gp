function dg_plotcomposites(IMFavg, IMFsem, commonstart, composite, varargin)
% Plotz the mean +/- 2 SEM and also the overlay of all composites,
% clickable for job index

% Options:
% 'period', period - sets sample period; overrides value in plotdata.
% 'plotdata', plotdata - uses <plotdata> returned by lfp_hht to plot event
%   markers and set sample period.
% 'ovr'- plots overlay of the IMFs from all the jobs (EMD param combos)

%$Rev: 75 $
%$Date: 2010-08-20 17:42:17 -0400 (Fri, 20 Aug 2010) $
%$Author: dgibson $


global lfp_SamplePeriod

ovrflag = false;
periodval = [];
plotdata = [];

argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'ovr'
            ovrflag = true;
        case 'period'
            argnum = argnum + 1;
            periodval = varargin{argnum};
        case 'plotdata'
            argnum = argnum + 1;
            plotdata = varargin{argnum};
        otherwise
            error('dg_plotcomposites:badoption', ...
                ['The option "' dg_thing2str(varargin{argnum}) '" is not recognized.'] );
    end
    argnum = argnum + 1;
end

if ~isempty(periodval)
    period = periodval;
elseif ~isempty(plotdata)
    period = plotdata.period;
else
    period = lfp_SamplePeriod;
end

timepts = period * (commonstart:size(IMFavg,2)+commonstart-1);

if ovrflag
    hF = figure;
    hA = axes('Parent', hF);
    hL = plot(hA, timepts, composite');
    set(hA, 'XGrid', 'on', 'YGrid', 'on');
    xlabel(hA, 'Time, s');
    for k = 1:size(composite,1)
        set(hL(k), 'ButtonDownFcn', sprintf('disp(''job %d'');', k));
    end
    if ~isempty(plotdata)
        lfp_plotEvtMarkers( hA, ...
            plotdata.trialevts(plotdata.evts2mark,:), ...
            'reftime', plotdata.reftime, ...
            'bigevts', plotdata.bigevts );
        lfp_createFigTitle(hA, [],[],[],[],[], 'figtitle', plotdata.figtitle);
    end
end

hF = figure;
hA = axes('Parent', hF, 'NextPlot', 'add');
plot(hA, timepts, ...
    [(IMFavg + 2*IMFsem)'  (IMFavg - 2*IMFsem)'], 'Color', [.7 .7 .7]);
plot(hA, timepts, IMFavg, 'k');
set(hA, 'XGrid', 'on', 'YGrid', 'on');
xlabel(hA, 'Time, s');
xlabel(hA, 'Time, s');
if ~isempty(plotdata)
    lfp_plotEvtMarkers( hA, ...
        plotdata.trialevts(plotdata.evts2mark,:), ...
        'reftime', plotdata.reftime, ...
        'bigevts', plotdata.bigevts );
    lfp_createFigTitle(hA, [],[],[],[],[], 'figtitle', plotdata.figtitle);
end

