function hL = lfp_plotEvtMarkers(hA, evtmatrix, varargin)
%hL = lfp_plotEvtMarkers(hA, evtmatrix, mode)
% Plot single or aggregated event marker(s) with detail strings.  Events
% that have NaN or 0 for the event ID are plotted with the default color
% and an empty event string.  (NaN should generally be reserved for events
% that do not exist in lfp_Events.)  Event markers are plotted as vertical
% lines that run the full height of the y axis.

%INPUTS
% hA: axes handle
% evtmatrix: a list of events in lfp_Events format
%
%OUTPUTS
% hL: a vector of line handles for the event markers
%
%OPTIONS
% 'bigevts', bigevtcodes - as in lfp_SpikeAnalysis. It is assumed that the
%   smallest ID in 'bigevts' is greater than the largest ID listed in
%   lfp_EventColors or lfp_EventNames.
% 'single' - plot individual event data for each row of <evtmatrix> (this
%   is the default)
% 'stats' - plot the "standard" set of 'evtavg' statistics, aggregated by
%   event IDs that actually appear in (evtmatrix(:,2)).  Does not handle
%   NaN event IDs well. Assumes that the timestamp of the reference event
%   in each trial has already been subtracted from the timestamps in
%   <evtmatrix>.
% 'reftime', reftime - plots events relative to <reftime>, but uses the
% 	absolute timestamps in the clickable detail strings.
%
%NOTES
% If lfp_EventDefaultColor is empty, assigns a value


%$Rev: 396 $
%$Date: 2019-02-15 18:13:19 -0500 (Fri, 15 Feb 2019) $
%$Author: dgibson $

global lfp_SelectedEventIDs lfp_SetupName %#ok<NUSED>

if length(hA) > 1
    error('lfp_plotEvtMarkers:nonscalar', ...
        '<hA> must be a scalar.');
end

bigevts = {};
bigevtIDs = [];
mode = 'single';
reftime = 0;
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'bigevts'
            argnum = argnum + 1;
            bigevts = varargin{argnum};
            try
                if isempty(bigevts)
                    bigevtIDs = [];
                else
                    bigevtIDs = cell2mat(bigevts(:,1));
                end
            catch %#ok<CTCH>
                error('lfp_plotEvtMarkers:badbigevts', ...
                    'Value for <bigevts> is badly formatted.' );
            end
        case {'single' 'stats'};
            mode = varargin{argnum};
        case 'reftime'
            argnum = argnum + 1;
            reftime = varargin{argnum};
        otherwise
            error('lfp_plotEvtMarkers:badoption', ...
                'The option "%s" is not recognized.', ...
                dg_thing2str(varargin{argnum}) );
    end
    argnum = argnum + 1;
end

if isempty(lfp_SelectedEventIDs)
    warning('lfp_plotEvtMarkers:setup', ...
        'Loading setup file');
    lfp_loadSetup;
end

set(hA, 'NextPlot', 'add');

switch mode
    case 'single'
        hL = NaN(size(evtmatrix,1), 1);
        for evtrow = 1:size(evtmatrix,1)
            evtid = evtmatrix(evtrow,2);
            evttime = evtmatrix(evtrow,1);
            if (evtid <= length(lfp_SelectedEventIDs)) ...
                    && lfp_SelectedEventIDs(evtid) ...
                    || ismember(evtid, bigevtIDs) ...
                    || isnan(evtid)
                % Either the event ID is selected, or it is a
                % specified big event ID.
                [evtname, evtcolor, evtcolorstr] = ...
                    lfp_getEvtProps(evtid, bigevts);
                hL(evtrow) = plot(hA, ...
                    [ evttime evttime ] - reftime, ...
                    get(hA, 'YLim'), 'Color', evtcolor );
                detailstr = sprintf( ...
                    '\\nTimestamp=%.6f\\nLineColor="%s"\\nEventName="%s"\\nEventID=%.0f (0x%X)', ...
                    evttime, evtcolorstr, evtname, ...
                    evtid, evtid );
                set(hL(evtrow), ...
                    'ButtonDownFcn', ...
                    ['fprintf(1,''' detailstr '\n'')'] );
            end
        end
    case 'stats'
        evtIDs = unique(evtmatrix(:,2))';
        hL = NaN(length(evtIDs), 1);
        for evtrow = 1:length(evtIDs)
            evtid = evtIDs(evtrow);
            [evtname, evtcolor, evtcolorstr] = lfp_getEvtProps(evtid, bigevts);
            TS = evtmatrix(evtmatrix(:,2)==evtid, 1);
            medianTS = median(TS);
            prct = prctile(TS, [5, 95, 25, 75]);
            detailstr = sprintf('\\nevtid=%d, "%s", color=%s approx stats\\nmedian=%6.4g, mean=%6.4g, SD=%6.4g\\nfirst quartile=%6.4g, last quartile=%6.4g\\n5th percentile=%6.4g, 95th percentile=%6.4g', ...
                evtid, evtname, evtcolorstr, medianTS, mean(TS), std(TS), prct(3), prct(4), prct(1), prct(2));
            hL(evtrow) = plot(hA, [ medianTS medianTS ], ...
                get(hA, 'YLim'), ...
                'Color', evtcolor );
            set(hL(evtrow), ...
                'ButtonDownFcn', ...
                ['fprintf(1,''' detailstr '\n'')'] );
        end
    otherwise
        error('lfp_plotEvtMarkers:badmode', ...
            'Unrecognized mode "%s"', dg_thing2str(mode));
end
end



