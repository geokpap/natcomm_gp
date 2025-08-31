function cells = dg_strsplit(filename)
% Returns contents of a tab-delimited text file formatted as a cell string
% array.
%INPUTS
% filename: path to a tab-delimited text file.
%OUTPUTS
% cells: cell string array with one row for each line in the file, and one
%   column for each tab-delimited column on the longest line in the file.
%NOTES
% This is incredibly stupid, but I had to replace 'strsplit' with my own
% version to this work correctly in R2022b ('strsplit' no longer returns
% empty cells when there are two delimiters in a row).

%$Rev: 306 $
%$Date: 2023-10-05 16:32:32 -0400 (Thu, 05 Oct 2023) $
%$Author: dgibson $

tab = sprintf('\t');
newline = sprintf('\n');
bytes = fileread(filename);
lines = strsplit(bytes, newline);
if isempty(lines{end})
    lines(end) = [];
end
cells = cell(length(lines), 0);
for linenum = 1:length(lines)
    values = strsplit(lines{linenum}, tab);
    cells(linenum, 1:length(values)) = values;
end
