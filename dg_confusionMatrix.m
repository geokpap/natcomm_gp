function [matrix, groupIDs] = dg_confusionMatrix(group, pred, varargin)
%[matrix, groupIDs] = dg_confusionMatrix(group, pred)
% Constructs the confusion matrix from results of a classifier (decoder)
% run.
%INPUTS
% group: the category that each trial actually belongs to.
% pred: the category that was predicted for each trial by the classifier.
%OUTPUTS
% matrix: number of trials that were actually in the group specified by the
%   row number and were classified in the group specified by the column.
% groupIDs: labels for the groups, sorted in the same order as the rows and
%   columns of <matrix>.
%OPTIONS
% 'order', order - <matrix> is ordered according to <groupIDs>, and
%   <groupIDs> is normally in alphabetic order.  This option allows you to
%   re-order <groupIDs> before constructing <matrix>; e.g. (...'order',
%   length(unique(group)):-1:1) would put groupIDs in reverse alphabetic
%   order.
%NOTES
% Use dg_plotConfusionMatrix to display the results.  <group> and <pred>
% are as in:
%   [hitrate, pred] = dg_leaveOneOut(data, group);

%$Rev: 186 $
%$Date: 2013-12-17 19:38:58 -0500 (Tue, 17 Dec 2013) $
%$Author: dgibson $

groupIDs = unique(group);
matrix = NaN(length(groupIDs));

argnum = 1;
while argnum <= length(varargin)
    if ischar(varargin{argnum})
        switch varargin{argnum}
            case 'order'
                argnum = argnum + 1;
                order = varargin{argnum};
                if ~isnumeric(order)
                    error('dg_confusionMatrix:order', ...
                        'The value for the ''order'' option must be numeric');
                end
                groupIDs = groupIDs(order);
            otherwise
                error('dg_confusionMatrix:badoption', ...
                    'The option %s is not recognized.', ...
                    dg_thing2str(varargin{argnum}));
        end
    else
        error('dg_confusionMatrix:badoption2', ...
            'The value %s occurs where an option name was expected', ...
            dg_thing2str(varargin{argnum}));
    end
    argnum = argnum + 1;
end


for realgrpidx = 1:length(groupIDs)
    isingroup = ismember(group, groupIDs(realgrpidx));
    for assignedgrpidx = 1:length(groupIDs)
        matrix(realgrpidx, assignedgrpidx) = ...
            sum(ismember(pred(isingroup), groupIDs(assignedgrpidx)));
    end
end
