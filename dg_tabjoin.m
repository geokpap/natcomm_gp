function joinedcells = dg_tabjoin(file1, file2, fieldnums1, fieldnums2, ...
    hasheader)
% Reads two tab-delimited text files and joins them like RDBMS tables on
% specified key fields.
%INPUTS
% file1, file2: paths to two tab-delimited text files.
% fieldnums1, fieldnums2: column numbers of a set of fields in each file
%   that define unique keys to use as the basis of the join.  The fields
%   need not be in the same order in both files, but they must be specified
%   in the same order in <fieldnums1>, <fieldnums2>.
% hasheader: a two element logical vector that specifies whether there is a
%   header for each file.
%OUTPUTS
% joinedcells: cell string array with one row for each line in the joined
%   table, and one column for each column in the joined table.  If either
%   input file has a header, then <joinedcells> also has a header.
%NOTES
%   The combinations of values in the key <fieldnums> must be unique across
% rows within each file.  If any combination of values in <fieldnums> in
% the second file does not exist in the first file, a warning is raised and
% the row in <file2> with the unmatched key fields is discarded.
%   The order of columns in the output is equal to the order of the
% columns in the first file plus the order of the columns excluding the key
% fields in the second file.  If the first file has no header, then the key
% fields in the output file will have empty headers.
%   Iteration is done over rows of the second file, so it will be faster to
% put the file with fewest rows second.
%   No assumption is made about the ordering of either input file.

%$Rev: 290 $
%$Date: 2022-04-01 21:48:18 -0400 (Fri, 01 Apr 2022) $
%$Author: dgibson $

if numel(fieldnums1) ~= numel(fieldnums2)
    error('dg_tabjoin:numfieldnums', ...
        '<fieldnums1>, <fieldnums2> must both have the same number of elements.');
end
numkeyflds = numel(fieldnums1);
header = {};
cells1 = dg_tabread3(file1);
if hasheader(1)
    header = cells1(1, :);
    cells1(1, :) = [];
end
keys1 = cell(size(cells1, 1), 1);
for rownum = 1:size(cells1, 1)
    keys1{rownum} = [cells1{rownum, fieldnums1}];
end
if size(unique(keys1), 1) ~= size(cells1, 1)
    error('dg_tabjoin:nonunique1', ...
        'The key fields in file1 are not unique.');
end
cells2 = dg_tabread3(file2);
fields2append = setdiff(1:size(cells2, 2), fieldnums2);
if hasheader(2)
    if isempty(header)
        header = ...
            [repmat({''}, 1, size(cells1, 2)) cells2(1, fields2append)];
    else
        header = [header cells2(1, fields2append)];
    end
    cells2(1, :) = [];
else
    if ~isempty(header)
        header = [header repmat({''}, 1, length(fields2append))];
    end
end
keys2 = cell(size(cells2, 1), 1);
for rownum = 1:size(cells2, 1)
    keys2{rownum} = [cells2{rownum, fieldnums2}];
end
if size(unique(keys2), 1) ~= size(cells2, 1)
    error('dg_tabjoin:nonunique2', ...
        'The key fields in file2 are not unique.');
end
joinedcells = cell( size(cells2, 1), ...
    size(cells1, 2) + size(cells2, 2) - numkeyflds );
for rownum = 1:size(cells2, 1)
    istarget = ismember(keys1, keys2{rownum});
    if ~any(istarget)
        warning('dg_tabjoin:nokey', ...
            'File1 does not contain key %s', keys2{rownum});
    else
        joinedcells(rownum, 1:size(cells1, 2)) = cells1(istarget, :);
        joinedcells(rownum, size(cells1, 2) + (1:length(fields2append))) ...
            = cells2(rownum, fields2append);
    end
end
% Remove any rows in <joinedcells> where the key fields could not be found:
joinedcells(cellfun(@isempty, joinedcells(:, fieldnums1(1))), :) = [];
% Prepend header, if any:
if ~isempty(header)
    joinedcells = [ header
        joinedcells ];
end
