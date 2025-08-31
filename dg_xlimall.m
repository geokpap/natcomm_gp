function dg_xlimall(hF, limits)
%dg_xlimall(hF, limits)
%  Sets XLim to <limits> on every axes in the figure whose handle is <hF>,
%  provided that its UserData property is empty (this is a quick and dirty
%  way to skip legend axes, and can be made more specific if necessary).
%  <limits> must be of the appropriate type and value to be assigned to the
%  'XLim' property of an axes.  For example, dg_xlimall(gcf, [0 50]) will
%  set all x-axes on the "current figure" to 0-50.

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

axeses = findobj(hF, 'Type', 'axes');
for k = 1:length(axeses)
    if isempty(get(axeses(k), 'UserData'))
        set(axeses(k), 'XLim', limits);
    end
end