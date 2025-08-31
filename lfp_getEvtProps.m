function [name, color, colorstr] = lfp_getEvtProps(evtid, bigevts)
%INPUTS
% evtid: index into lfp_EventColors, etc.
% bigevts: value suitable for <bigevtcodes> as in lfp_SpikeAnalysis or 
%   lfp_plotEvtMarkers 'bigevts' option.

%OUTPUTS
% name: string suitable for embedding in click string.
% color: suitable for use as the value of the 'Color' option to 'plot'.
% colorstr: string suitable for embedding in click string.

%$Rev: 396 $
%$Date: 2019-02-15 18:13:19 -0500 (Fri, 15 Feb 2019) $
%$Author: dgibson $

global lfp_EventDefaultColor lfp_EventColors lfp_EventNames
if isempty(lfp_EventDefaultColor)
    lfp_EventDefaultColor = 'b';
end
if evtid == 0 || isnan(evtid)
    color = lfp_EventDefaultColor;
    name = '';
else
    color = '';
    if evtid <= length(lfp_EventColors)
        color = lfp_EventColors{evtid};
    elseif ~isempty(bigevts) && ismember(evtid, cell2mat(bigevts(:,1)))
        for k = 1:size(bigevts,1)
            if ismember(evtid, bigevts{k,1})
                name = sprintf('bigevtid=%.0f', evtid); 
                color = bigevts{k,2};
                if isempty(color)
                    color = lfp_EventDefaultColor;
                end
                break
            end
        end
    end
    if evtid <= length(lfp_EventNames)
        name = lfp_EventNames{evtid};
    elseif isempty(bigevts) || ~ismember(evtid, cell2mat(bigevts(:,1)))
        name = '';
    end
end
if isempty(color)
    color = lfp_EventDefaultColor;
end
if ischar(color)
    colorstr = color;
else
    colorstr = mat2str(color);
end
end