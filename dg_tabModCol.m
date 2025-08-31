function cells = dg_tabModCol(funch, cells, colnum, hasheader)
% Creates a modified copy of a column in a 'dg_tabread3' cell array.
%INPUTS
% funch: function handle to apply to each individual cell in <colnum> to
%   perform the modifications.
% cells: as returned by 'dg_tabread3'.
% colnum: column number to copy and modify.  The modified copy is appended
%   at the end of <cells>.
% hasheader: a two element logical vector that specifies whether there is a
%   header for each file.
%OUTPUTS
% cells: modified copy of input argument <cells>.
%NOTES
% If there is a header, cells{1, colnum} is left unmodified, and its value
% with 'mod' appended is used as the header of the new column.

%$Rev: 292 $
%$Date: 2022-04-28 17:28:01 -0400 (Thu, 28 Apr 2022) $
%$Author: dgibson $

lastcol = size(cells, 2);
if hasheader(1)
    startrow = 2;
    cells{1, lastcol+1} = [cells{1, colnum} 'mod'];
else
    startrow = 1;
end
for rownum = startrow:size(cells, 1)
    cells{rownum, lastcol+1} = feval(funch, cells{rownum, colnum});
end

