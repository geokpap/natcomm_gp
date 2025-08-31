function cells = dg_tabread3(filename, varargin)
% Returns contents of a tab-delimited text file formatted as a cell string
% array.
%INPUTS
% filename: path to a tab-delimited text file.
%OUTPUTS
% cells: cell string array with one row for each line in the file, and one
%   column for each tab-delimited column on the longest line in the file.
%OPTIONS
% 'delim', delim - <delim> replaces tab in the line parsing.  It can be any
%   single character.
%NOTES
% Returns an empty cell (i.e., {''}) when there are two delimiters in a
% row.
%   Note that DOS-formatted text files will include a <CR> at the end of
% the last cell in each row, and will thus produce another DOS-formatted
% text file if <cells> is converted back to text using Unix formatting.

%$Rev: 306 $
%$Date: 2023-10-05 16:32:32 -0400 (Thu, 05 Oct 2023) $
%$Author: dgibson $

delim = sprintf('\t');
newline = sprintf('\n');

argnum = 0;
while true
    argnum = argnum + 1;
    if argnum > length(varargin)
        break
    end
    if ~ischar(varargin{argnum})
        continue
    end
    switch varargin{argnum}
        case 'delim'
            argnum = argnum + 1;
            delim = varargin{argnum};
        otherwise
            error('dg_tabread3:badoption', ...
                'The option %s is not recognized.', ...
                dg_thing2str(varargin{argnum}));
    end
end
if numel(delim) > 1
    error('dg_tabread3:delim', ...
        '<delim> must be a single character.');
end

bytes = fileread(filename);
isCR = bytes == 13;
if any(isCR)
    warning('dg_tabread3:CR', 'Removing <CR> characters.');
    bytes(isCR) = [];
end
lines = strsplit(bytes, newline, 'CollapseDelimiters', false);
if isempty(lines{end})
    lines(end) = [];
end
cells = cell(length(lines), 0);
for linenum = 1:length(lines)
    values = strsplit(lines{linenum}, delim, 'CollapseDelimiters', false);
    cells(linenum, 1:length(values)) = values;
end
