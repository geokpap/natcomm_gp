function hL = dg_plotContours(C, hA, mywidth, mycolor)
%dg_plotContours(C)
%dg_plotContours(C, hA)
%dg_plotContours(C, hA, width)
%dg_plotContours(C, hA, width, color)
% INPUTS
%   C: contour data as returned by Matlab func 'contourc'
%   hA: an axes into which to plot the contours; if given they are plotted
%       as 2 point white lines, and if not, then a new figure is created
%       and they are plotted as 0.5 point black lines.
%   width: overrides the default widths described above.
%   color: overrides the default colors described above.
%OUTPUT
% hL: row vector of handles to the line objects plotted.

%$Rev: 221 $
%$Date: 2015-06-18 16:53:09 -0400 (Thu, 18 Jun 2015) $
%$Author: dgibson $
if nargin < 2
    figure;
    hA = axes('NextPlot', 'add');
    if nargin < 4
        mycolor = 'k';
    end
    if nargin < 3
        mywidth = 0.5;
    end
else
    set(hA, 'NextPlot', 'add');
    if nargin < 4
        mycolor = 'w';
    end
    if nargin < 3
        mywidth = 2.0;
    end
end
idx = 1;
hL = zeros(1, size(C,2));
while idx <= size(C,2)
    hL(idx) = plot(hA, C(1, idx+1:idx+C(2,idx)), C(2, idx+1:idx+C(2,idx)), ...
        'LineWidth', mywidth, ...
        'Color', mycolor);
    idx = idx + C(2,idx) + 1;
end