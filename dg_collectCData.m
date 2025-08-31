function [cdata, xvals, yvals] = dg_collectCData(hF)
%[cdata, xvals, yvals] = dg_collectCData(hF)
% Concatenates the CData from all the panels in the figure <hF>, which is
% typically a pasteup from lfp_makepasteup.  It is assumed that the images
% will be ordered from right to left in the list of images contained in
% <hF>, and if any of them is a colorbar it is ignored.  All panels must
% have the same y values.

%$Rev: 108 $
%$Date: 2011-05-10 01:04:31 -0400 (Tue, 10 May 2011) $
%$Author: dgibson $

hI = findobj(hF, 'Type', 'image');
cdata = {};
xvals = {};
yvals = [];
for imgnum = 1:length(hI)
    if isequal(get(get(hI(imgnum), 'Parent'), 'Tag'), 'Colorbar')
        % just skip it
    else
        cdata{1,end+1} = get(hI(imgnum), 'CData');
        xvals{1,end+1} = get(hI(imgnum), 'XData');
        ytemp = get(hI(imgnum), 'YData');
        if ~isempty(yvals) && ~isequal(yvals, ytemp)
            error('dg_collectCData:badyvals2', ...
                'Images number %d haa different yvals from preceding', ...
                imgnum);
        else
            yvals = ytemp;
        end
    end
end
cdata = cell2mat(cdata(end:-1:1));
xvals = cell2mat(xvals(end:-1:1));

