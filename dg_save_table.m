function dg_save_table(header, keys, data, FileName)
%dg_save_table saves a table as a tab-delimited spreadsheet.
%  Open <FileName>.xls for writing and save the table <data> with <keys>
%  inserted as the first column.  <header> is a cell vector of
%  mat2str-compatible elements for the header row (that includes both
%  strings and doubles).  <FileName> will be silently overwritten if it
%  already exists.  If <FileName> is not given or is empty, then a GUI is
%  presented for specifying the FileName.
%
%  Obviously, <keys> must have the same number of elements as <data> has
%  rows, and <header> must have the same number of elements as <data> has
%  columns. Element (1,1) of the output table is left empty.  If <header>
%  is [], then no header row is produced.  If <keys> is [], then the first
%  column is empty.
%
%  <keys> and <data> are both assumed to be numeric.

%$Rev: 153 $
%$Date: 2012-07-17 18:40:53 -0400 (Tue, 17 Jul 2012) $
%$Author: dgibson $

if nargin < 4
    FileName = [];
end
if ~isempty(header) && ~((size(header,1) == 1) || (size(header,2) == 1))
    error('dg_save_table:badheader2', ...
        '<header> must be vector' );
end
if ~isempty(keys) && ~((size(keys,1) == 1) || (size(keys,2) == 1))
    error('dg_save_table:badkeys2', ...
        '<keys> must be vector' );
end
if ~isempty(header) && (size(data,2) ~= length(header))
    error('dg_save_table:badheader', ...
        '<header> must have same number of elements as <data> has columns');
end
if ~isempty(keys) && (length(keys) ~= size(data,1))
    error('dg_save_table:badkeys', ...
        '<keys> must have same number of elements as <data> has rows');
end
if isempty(FileName)
    [FileName,PathName] = uiputfile('spreadsheet.xls');
else
    [PathName,name] = fileparts(FileName);
    FileName = [name '.xls'];
end
if isequal(FileName, 0) return
end
fid = fopen(fullfile(PathName, FileName), 'w');
if fid == -1
    error('dg_save_table:badfilename', ...
        'Could not open file %s for writing', fullfile(FileName,PathName) );
end
if ~isempty(header)
    for col = 1:length(header)
        fprintf(fid, '\t%s', mat2str(header{col}));
    end
    fprintf(fid, '\n');
end
for row = 1:size(data,1)
    if ~isempty(keys)
        fprintf(fid, '%d', keys(row));
    end
    fprintf(fid, '\t%d', data(row,:));
    fprintf(fid, '\n');
end
fclose(fid);