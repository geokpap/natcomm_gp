function num = dg_spaceread(filename)
%DG_SPACEREAD reads space-delimited text spreadsheets
%num = dg_spaceread(filename)
% Behaves similarly to num = xlsread(filename), except that <filename> must
% be the pathname of a space-delimited text file:  dg_spaceread ignores
% leading rows or columns of text; however, if a cell not in a leading row
% or column contains text, dg_spaceread puts a NaN in its place in the
% return array, num.  Tolerates unequal numbers of fields per line. Note
% that because any number of delims in a row is considered to be one
% delimiter, it is impossible for a cell to be empty in the text file.
% Ignores any lines that begin with '#' (the "pound" character).

% Tested only on a file with spaces at beginnings of lines but not ends.

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

rowallocsize = 8192;
colallocsize = 5;
num = zeros(rowallocsize); % pre-allocate for speed
linenum = 0;
numcols = 0;

fid = fopen(filename);
if fid == -1
    error('Could not open %s', filename);
end

line = fgetl(fid);
while ~isequal(line, -1)
    if line(1)~='#'
        linenum = linenum + 1;
        delims = regexp(line, ' +');
        starts = regexp(line, '[^ ]+');
        % Remove leading delimiter if there is one:
        if delims(1) == 1
            delims(1) = [];
        end
        % Add trailing delimiter if there is none:
        if delims(end) <= starts(end)
            delims(end+1) = length(line) + 1;
        end
        numcols = max(numcols, length(delims));
        % Allocate more storage if needed:
        if linenum > size(num,1)
            num = [ num; zeros(rowallocsize, size(num,2)) ];
        end
        if numcols > size(num,2)
            num = [ num zeros(size(num,1), colallocsize) ];
        end
        % For each field, there is now one start and one delim, and
        % starts(k) < delims(k).
        % convert text to number:
        for delimnum = 1:length(delims)
            value = str2num(line(starts(delimnum) : delims(delimnum)-1));
            if isempty(value)
                num(linenum, delimnum) = NaN;
            else
                num(linenum, delimnum) = value;
            end
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

fclose(fid); 