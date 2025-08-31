function [hI, hCB] = dg_plotConfusionMatrix(hDisp, matrix, groupIDs)
%hF = dg_plotConfusionMatrix(matrix, groupIDs)
% Creates a pseudocolor plot of a confusion matrix.  The matrix is
% always normalized so that each row totals to 1 before plotting.
%INPUTS
% hDisp: graphics handle to the object in which the pseudocolor plot will
%   be created.  If empty ([] or ''), then a new figure is created.
% matrix: the confusion matrix, in Predicted X Actual format (as returned
%   by dg_confusionMatrix).
% groupIDs: vector of identifiers to be used as tick labels on the axes of
%   the pseudocolor plot (as returned by dg_confusionMatrix).
%OUTPUTS
% hI: axes handle to newly created pseudocolor image.
% hCB: axes handle to newly created colorbar that goes with <hI>.
%NOTES
% <matrix> and <groupIDs> are as returned by:
%   [matrix, groupIDs] = dg_confusionMatrix(group, class);

%$Rev: 179 $
%$Date: 2013-09-06 13:09:46 -0400 (Fri, 06 Sep 2013) $
%$Author: dgibson $

if isempty(hDisp)
    hDisp = figure;
end

switch get(hDisp, 'Type')
    case 'axes'
        hA = hDisp;
    case 'figure'
        hA = axes('Parent', hDisp);
    otherwise
        error('dg_plotConfusionMatrix:unknowndisp', ...
            'Unknown display object type: %s', get(hDisp, 'Type'));
end

normmatrix = matrix ./ repmat(sum(matrix,2), 1, size(matrix,2));
hI = imagesc(normmatrix, 'Parent', hA);
axis(hA, 'equal');
caxis(hA, [0 1]);
hCB = colorbar('peer', hA);
set(get(hCB,'YLabel'), 'String', 'Conditional probability of prediction');
ylabel(hA, 'Fraction of Trials Predicted');
xlabel(hA, 'Predicted Trial Type');
ylabel(hA, 'Actual Trial Type');
set(hA, 'XLim', [0.5 length(groupIDs)+0.5])
set(hA, 'YLim', [0.5 length(groupIDs)+0.5])
set(hA, 'XTick', 1:length(groupIDs));
set(hA, 'YTick', 1:length(groupIDs));
set(hA, 'XTickLabel', groupIDs);
set(hA, 'YTickLabel', groupIDs);
