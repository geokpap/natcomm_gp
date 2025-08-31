function [num, rowhead, colhead] = dg_tabreadfix(filename)
%DG_TABREADFIX reads tab-delimited text spreadsheets in fixed-point format.
%[num, txt] = dg_tabreadfix(filename)
% Identical to dg_tabread except that sscanf is used instead of str2num;
% this means that strings starting with numbers do convert to numbers
% rather than to NaN.  Also, if the return values <rowhead> and <colhead>
% are used, then they contain cell string arrays of the character data in
% the first column and first row respectively.  Also, the portion of code
% during which the file is open is wrapped in a try...catch that closes the
% file.
%
% Behaves similarly to num = xlsread(filename), except that <filename> must
% be the pathname of a tab-delimited text file:  dg_tabreadfix ignores leading
% rows or columns of text; however, if a cell not in a leading row or
% column is empty or contains text, dg_tabreadfix puts a NaN in its place in
% the return array, num.  Tolerates unequal numbers of tabs per line.

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

allocsize = 128;
num = zeros(allocsize); % pre-allocate for speed
if nargout > 1
    rowhead = cell(allocsize,1); % pre-allocate for speed
    if nargout > 2
        colhead = {};
    end
end
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
        tabs = regexp(line, '\t');
        tabs(end+1) = length(line) + 1;
        numcols = max(numcols, length(tabs));
        if linenum == 1 && nargout > 2
            colhead{1} = line(1 : tabs(1)-1);
            for tabnum = 2:length(tabs)
                colhead{tabnum} = line(tabs(tabnum-1)+1 : tabs(tabnum)-1);
            end
        end
        % Allocate more storage if needed:
        if linenum > size(num,1)
            num = [ num; zeros(allocsize, size(num,2)) ];
            if nargout > 1
                rowhead = [ rowhead; cell(allocsize,1) ];
            end
        end
        if numcols > size(num,2)
            num = [ num zeros(size(num,1), allocsize) ];
        end
        % convert text to number:
        if nargout > 1
            rowhead{linenum} = line(1 : tabs(1)-1);
        end
        value = sscanf(line(1 : tabs(1)-1), '%f');
        if isempty(value)
            num(linenum, 1) = NaN;
        else
            num(linenum, 1) = value;
        end
        for tabnum = 2:length(tabs)
            value = sscanf(line(tabs(tabnum-1)+1 : tabs(tabnum)-1), '%f');
            if isempty(value)
                num(linenum, tabnum) = NaN;
            else
                num(linenum, tabnum) = value;
            end
        end
        line = fgetl(fid);
    end
    fclose(fid);
catch
    msg = lasterror;
    fclose(fid);
    error(msg);
end
% Trim off unused allocated storage:
num = num(1:linenum, 1:numcols);
if nargout > 1
    rowhead = rowhead(1:linenum);
end
% Trim off empty leading rows and columns:
while ~isempty(num) && all(isnan(num(1,:)))
    num(1,:) = [];
end
while ~isempty(num) && all(isnan(num(:,1)))
    num(:,1) = [];
end
