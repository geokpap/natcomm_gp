function A = dg_tolSetdiff(A, B, tol)
%result = dg_tolSetdiff(A, B, tol)
%   Returns the same result as setdiff(A,B,'rows') except that values of
%   A(:,end) and B(:,end) that differ by less than tol are considered to
%   be equal.  Note that if tol = 0, then no elements are eliminated.

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

A = unique(A, 'rows');
B = unique(B, 'rows');
if size(A,2) > 1
    samestartcols = ismember(A(:, 1:end-1), B(:, 1:end-1), 'rows');
else
    samestartcols = true(size(A,1), 1);
end
rows2delete = [];
ismatch = false(size(A,1), 1);
for rowA = find(samestartcols)'
    if size(B,2) > 1
        istargetB = ismember(B(:, 1:end-1), A(rowA, 1:end-1), 'rows');
    else
        istargetB = true(size(B,1), 1);
    end
    rowBlist = find(istargetB);
    ismatch(rowA) = any(abs(B(rowBlist, end) - A(rowA, end)) < tol);
end
A(ismatch, :) = [];
