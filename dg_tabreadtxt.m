function txt = dg_tabreadtxt(filename)
%txt = dg_tabreadtxt(filename)
% Identical to dg_tabread except that all fields are treated as literal
% text. Behaves similarly to [~, txt] = xlsread(filename), except that
% <filename> must be the pathname of a tab-delimited text file, and all
% fields are treated as literal text.

%$Rev:  $
%$Date:  $
%$Author: dgibson $

txt = {};
linenum = 0;
numcols = 0;

fid = fopen(filename);
if fid == -1
    error('Could not open %s', filename);
end

line = fgetl(fid);
while ~isequal(line, -1)
    linenum = linenum + 1;
    tabs = [regexp(line, '\t') length(line) + 1];
    numcols = max(numcols, length(tabs));
    txt{linenum, 1} = line(1 : tabs(1)-1); %#ok<AGROW>
    for tabnum = 2:length(tabs)
        txt{linenum, tabnum} = line(tabs(tabnum-1)+1 : tabs(tabnum)-1); %#ok<AGROW>
    end
    line = fgetl(fid);
end

fclose(fid);