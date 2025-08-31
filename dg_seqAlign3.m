function [sim, pairs] = dg_seqAlign3(v, w)
%DEG_SEQALIGN2 dynamic programming algorithm for LCS sequence alignment
%[sim, pairs] = dg_seqAlign3(v, w)

%[sim, pairs] = dg_seqAlign3(v, w)
% Modifed for speed from dg_seqAlign2 by DG 1/4/2007.
% <v> and <w> are arrays of elements that can be compared using '=='.
% <sim> is the accumulated similarity score between <v> and <w>.  <pairs>
% is a two-column array containing the indices of matched elements in <v>
% and <w>, respectively. 
%
% The values in <indel> are customarily negative, as are the off-diagonal
% elements of <match>.  The diagonal of <match> is customarily 1.  There is
% no special significance to the sign of the scores, however; what really
% matters is the difference between the maximum possible score (deltamax)
% and the score actually assigned to the edit operation.  To compute the
% Longest Common Subsequence, set the off-diagonal elements of <match> =
% -Inf and set indel(:) = 0.

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

v = reshape(v,[],1);
w = reshape(w,1,[]);

% Initialize arrays
s = uint16(zeros(numel(v)+1, numel(w)+1));  % cumulative similarities
p = int8(zeros(numel(v)+1, numel(w)+1));  % ptrs to predecessors: 1=del, 2=ins, 3=match

% Fill in non-boundaries
for i = 2:numel(v)+1
    for j = 2:numel(w)+1
        if v(i-1) == w(j-1)
            match = 1;
        else
            match = -Inf;
        end
        [s(i,j), p(i,j)] = max(...
            [ s(i-1,j), ...
                s(i,j-1), ...
                s(i-1,j-1) + match ]);
    end
end

sim = s(end,end);

pairs = zeros(max(numel(v), numel(w)), 2);
pairnum = size(pairs, 1);
while (i > 1 && j > 1)
    if p(i,j) == 3
        pairs(pairnum, :) = [i-1, j-1];
        pairnum = pairnum - 1;
        i = i - 1;
        j = j - 1;
    elseif p(i,j) == 2
        j = j - 1;
    elseif p(i,j) == 1
        i = i - 1;
    else
        error('dg_seqAlign3:badp', ...
            'Internal error.' );
    end
end
if pairnum
    pairs(1:pairnum, :) = [];
end

