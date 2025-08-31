function num = dg_tabread(filename, varargin)
%DG_TABREAD reads tab-delimited text spreadsheets
%INPUTS
% filename: absolute or relative path to file to read.
%OUTPUTS
% num: Behaves similarly to num = xlsread(filename), except that <filename>
%   must be the pathname of a tab-delimited text file:  dg_tabread ignores
%   leading rows or columns of text; however, if a cell not in a leading
%   row or column is empty or contains text, dg_tabread puts a NaN in its
%   place in the return array, num.  Tolerates unequal numbers of tabs per
%   line. "Text" and "numeric" strings are distinguished by passing the
%   string to str2double, and if the result is NaN, then the string was
%   text.
%OPTIONS
% 'delim', delim - <delim> replaces tab in the line parsing.  It can be any
%   single character.
% 'allocsize', allocsize - overrides default value of 512 for <allocsize>.
%   However, increasing <allocsize> may actually *increase* the running
%   time, if it is too large.
%NOTES
% The presence of even a single numerical field in a row or column causes
% that raw or column to be included in <num>, i.e. NOT to be treated as a
% header.  Also, if there is a header row containing fewer headers than
% there are columns of numerical data in the longest row in the file, then
% the missing columns in the header row are filled in with the value zero,
% and the entire row is therefore included in <num>, not trimmed off.
%   This function is very slow.  See 'dg_tabread2' and 'dg_tabread3'.

%$Rev: 306 $
%$Date: 2023-10-05 16:32:32 -0400 (Thu, 05 Oct 2023) $
%$Author: dgibson $

allocsize = 512;
delim = sprintf('\t');

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
        case 'allocsize'
            argnum = argnum + 1;
            allocsize = varargin{argnum};
        case 'delim'
            argnum = argnum + 1;
            delim = varargin{argnum};
        otherwise
            error('dg_tabread:badoption', ...
                'The option %s is not recognized.', ...
                dg_thing2str(varargin{argnum}));
    end
end
if numel(delim) > 1
    error('dg_tabread:delim', ...
        '<delim> must be a single character.');
end

num = zeros(allocsize); % pre-allocate for speed
linenum = 0;
numcols = 0;

fid = fopen(filename);
if fid == -1
    error('Could not open %s', filename);
end

try
    line = fgetl(fid);
    while ~isequal(line, -1)
        linenum = linenum + 1;
        tabs = strfind(line, delim);
        tabs(end+1) = length(line) + 1; %#ok<AGROW>
        numcols = max(numcols, length(tabs));
        % Allocate more storage if needed:
        if linenum > size(num,1)
            num = [ num; zeros(allocsize, size(num,2)) ]; %#ok<AGROW>
        end
        if numcols > size(num,2)
            num = [ num zeros(size(num,1), allocsize) ]; %#ok<AGROW>
        end
        % convert text to number:
        value = str2double(line(1 : tabs(1)-1));
        if isnan(value)
            num(linenum, 1) = NaN;
        else
            num(linenum, 1) = value;
        end
        for tabnum = 2:length(tabs)
            value = str2double(line(tabs(tabnum-1)+1 : tabs(tabnum)-1));
            if isnan(value)
                num(linenum, tabnum) = NaN;
            else
                num(linenum, tabnum) = value;
            end
        end
        line = fgetl(fid);
    end
    % Trim off unused allocated storage:
    num = num(1:linenum, 1:numcols);
    % Trim off empty leading rows and columns:
    while ~isempty(num) && all(isnan(num(1,:)))
        num(1,:) = [];
    end
    while ~isempty(num) && all(isnan(num(:,1)))
        num(:,1) = [];
    end
catch e
    fclose(fid);
    rethrow(e);
end

fclose(fid);