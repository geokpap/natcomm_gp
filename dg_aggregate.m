function [nums, cells] = dg_aggregate(nums, cells, cols)
%DG_AGGREGATE a cell array according to the contents of a descriptive
%numeric array.
%[nums, cells] = dg_aggregate(nums, cells, cols)

% <nums> is a numerical array, <cells> is a cell array with the same number
% of rows as <nums>.  Finds rows of <nums> that are identical in all
% columns *except* those listed in <cols>, removes them, and replaces them
% with a single row having NaN in columns <cols>.  Does a parallel removal
% of the same rows from <cells> that were removed from <nums>, and replaces
% them with a single row, each of whose elements is the vertical
% concatenation of the cell contents from the same column of the removed
% rows. This implies that the contents of the cells must be vertically
% concatenable.

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

includecols = setdiff(1:size(nums,2), cols);
agg_groups = cell(size(nums,1),1);  % list of lists of rows to aggregate
agg_cells = {};    % array of aggregated cell rows
agg_nums = [];     % array of aggregated numerical descriptors

% construct the aggregates
for row = 1 : (size(nums,1) - 1)
    % If <row> is already in an agg_group, skip it
    if ~any(cell2mat(dg_mapfunc( @ismember, ...
            repmat({row}, size(agg_groups)), agg_groups )))
        % look for a new agg_group containing <row>
        for row2 = (row+1) : size(nums,1)
            if isequalwithequalnans(...
                    nums(row, includecols), nums(row2, includecols) )
                if isempty(agg_groups{row})
                    agg_groups{row} = row;
                end
                agg_groups{row}(end+1,1) = row2;
            end
        end
        % if we found a new agg_group, add the aggregate to the list
        if ~isempty(agg_groups{row})
            agg_cells(end+1,:) = cells(agg_groups{row}(1),:);
            for k = 2:length(agg_groups{row})
                agg_cells(end,:) = ...
                    dg_mapfunc(@vertcat, ...
                    agg_cells(end,:), ...
                    cells(agg_groups{row}(k),:) );
            end
            agg_nums(end+1,:) = nums(row,:);
            agg_nums(end, cols) = NaN;
        end
    end
end

% remove the original rows that were agg_cells and append the aggregates
cells(cell2mat(agg_groups),:) = [];
cells = [cells; agg_cells];
nums(cell2mat(agg_groups),:) = [];
nums = [nums; agg_nums];

        
