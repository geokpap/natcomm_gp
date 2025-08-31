function dg_seqdisp(s1)
%dg_seqdisp(s1)
%   <s1> is a cell vector each of whose elements is a non-cell vector.
%   Regardless of the orientation of s1 or its elements, displays each cell
%   as a row vector.

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

for k = 1:length(s1)
    disp(reshape(s1{k}, 1, []));
end