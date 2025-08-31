function dg_ylimall(hF, limits)
%dg_ylimall(hF, limits)
%  Same as dg_xlimall but for the y axis.

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

axeses = findobj(hF, 'Type', 'axes');
for k = 1:length(axeses)
    if isempty(get(axeses(k), 'UserData'))
        set(axeses(k), 'YLim', limits);
    end
end