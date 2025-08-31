function binidx = dg_bindata(data, edges)
%INPUTS
% data: a 2-column array with x values in left column and y values in
%   right column.
% edges: defines bins similarly to Matlab function 'histc': a row of <data>
%   belongs to bin <k> if edges(k) <= data(i,1) < edges(k+1).  The
%   difference is that the last bin counts any values of x that are GREATER
%   THAN OR equal to edges(end), rather than the generally useless empty
%   list of values that EXACTLY match edges(end).
%OUTPUTS
% binidx: a matrix containing the same number of rows as <data> and
%   length(edges) columns.  Column <k> is a logical index into the rows of
%   <data> which is true for rows that belong to bin <k>.  In other words,
%   the set of all data points that belong to bin <k> is given by
%   data(binidx(:,k), :).

%$Rev: 123 $
%$Date: 2011-06-18 19:49:26 -0400 (Sat, 18 Jun 2011) $
%$Author: dgibson $

binidx = false(size(data,1), length(edges));
for k = 1 : (length(edges) - 1)
    binidx(:,k) = edges(k) <= data(:,1) & data(:,1) < edges(k+1);
end
binidx(:,end) = edges(end) <= data(:,1);
