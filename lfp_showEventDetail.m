function lfp_showEventDetail(eventindex)
%DEPRECATED

%lfp_showEventDetail(eventindex)
%  Displays information about the event on row <eventindex> in lfp_Events.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

warning('lfp_showEventDetail:deprecated', ...
    'This function is deprecated.')
lfp_declareGlobals;
if (lfp_Events(eventindex,2) == 0) || ...
        (lfp_Events(eventindex,2) > length(lfp_EventNames))
    linecolor = lfp_EventDefaultColor;
    eventname = '';
else
    eventname = lfp_EventNames{lfp_Events(eventindex,2)};
    if isempty(eventname)
        eventname = '';
    end
    linecolor = lfp_EventColors{lfp_Events(eventindex,2)};
    if isempty(linecolor)
        linecolor = '';
    end
end
fprintf(1, ...
    '\nTimestamp=%.6f\nLineColor="%s"\nEventName="%s"\nEventID=%d (0x%X)\n', ...
    lfp_Events(eventindex,1), ...
    mat2str(linecolor), eventname, ...
    lfp_Events(eventindex,2), lfp_Events(eventindex,2) ...
    );