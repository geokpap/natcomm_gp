function result = dg_tolIsmember(A, B, tol)
%result = dg_tolIsmember(A, B, tol)
%   Returns the same result as ismember(A,B) where A and B are both numeric
%   vectors, except that values of A and B that differ by strictly less
%   than tol are considered to be equal.  Note that if tol = 0, then
%   <result> is all false. It is assumed that A and B are already sorted in
%   ascending
%   order.
%NOTES
% The Matlab 'ismember' says it uses a method which sorts list, then
% performs binary search.  However, it calls a C function to do the actual
% binary search, so this won't be as fast.  However, we *can* take
% advantage of the fact that *both* lists are sorted.


%$Rev: 239 $
%$Date: 2016-03-10 19:30:16 -0500 (Thu, 10 Mar 2016) $
%$Author: dgibson $

result = false(size(A));
for aidx = 1:length(A)
    if aidx == 1
        bidx = dg_binsearch(B, A(aidx));
    else
        % since both lists are sorted, we kxnow that <bidx> has to be at
        % least as large as it was last iteration.
        bidx = dg_binsearch(B(bidx:end), A(aidx)) + bidx - 1;
    end
    if isempty(bidx)
        result(aidx) = A(aidx) < B(end) + tol;
    else
        result(aidx) = bidx <= length(B) && ...
            abs(B(bidx) - A(aidx)) < tol || ...
            (bidx > 1 && abs(B(bidx - 1) - A(aidx)) < tol);
    end
end
