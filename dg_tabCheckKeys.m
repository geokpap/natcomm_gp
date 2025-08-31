function [isunique, repeats] = dg_tabCheckKeys(cells, fieldnums)
% Verifies that the combination of key fields is unique across all records
% in a 'dg_tabread3' cell array.
%INPUTS
% cells: as returned by 'dg_tabread3'.
% fieldnums: column numbers of a set of fields that conjointly define
%   a key value that is unique for each row in the table.
%OUTPUTS
% isunique: <true> if the result of joining the key fields is unique for
%   each row in the file.
% repeats: if <~isunique>, then <repeats> is a two-column array containing
%   the row numbers of each repeated key and the row on which it is first
%   repeated.
%NOTES
% It is assumed that the header row, if any, has already been removed.

%$Rev: 292 $
%$Date: 2022-04-28 17:28:01 -0400 (Thu, 28 Apr 2022) $
%$Author: dgibson $

repeats = [];
keys1 = cell(size(cells, 1), 1);
for rownum = 1:size(cells, 1)
    keys1{rownum} = [cells{rownum, fieldnums}];
end
isunique = size(unique(keys1), 1) == size(cells, 1);
if ~isunique && nargout > 1
    for rownum = 1:size(cells, 1)
        isrepeat = ismember(keys1(rownum+1:end), keys1{rownum});
        if any(isrepeat)
            repeats(end+1, 1:2) = [rownum, find(isrepeat, 1) + rownum]; %#ok<AGROW>
        end
    end
end

