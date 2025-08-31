function dg_cell2tab(filename, C)
%dg_cell2tab(filename, C)
% Writes a tab-delimited text spreadsheet file containing the cell string
% array <C>.

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

fid = fopen(filename, 'w');
if isequal(fid, -1)
    error('dg_cell2tab:badfn', ...
        'Could not write to "%s"', filename);
end
try
    for row = 1:size(C,1)
        for col = 1 : (size(C,2) - 1)
            fprintf(fid, '%s\t', C{row,col});
        end
        fprintf(fid, '%s\n', C{row,end});
    end
catch
    fclose(fid);
    rethrow(lasterror);
end
fclose(fid);

            