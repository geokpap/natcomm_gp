function dg_alignPairedSubplots(k)
% operates on current figure
% for subplots produced using dg_subplot(... 'paired')
% k=1 aligns the upper of each pair to the lower.

%$Rev: 24 $
%$Date: 2009-03-31 21:51:08 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

if ~ismember(k, [1 2])
    error('<k> must be 1 or 2.');
end
appdata = getappdata(gcf);
appdatanames = fieldnames(appdata);
if ismember('SubplotGrid', appdatanames)
    grid = appdata.SubplotGrid;
else
    error('The current figure must contain subplots.');
end
for row = 1:2:size(grid,1)
    if row < size(grid,1)
        for col = 1:size(grid,2)
            if k == 1
                hA = dg_subplot(size(grid,1), size(grid,2), [row col]);
                hRef = dg_subplot(size(grid,1), size(grid,2), [row+1 col]);
            else
                hA = dg_subplot(size(grid,1), size(grid,2), [row+1 col]);
                hRef = dg_subplot(size(grid,1), size(grid,2), [row col]);
            end
            set(hA, 'XLim', get(hRef, 'XLim'));
        end
    end
end
