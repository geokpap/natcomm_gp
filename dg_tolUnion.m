function result = dg_tolUnion(A, B, tol)
%result = dg_tolUnion(A, B, tol)
%   Returns the same result as union(A,B,'rows') except that values of
%   A(:,end) and B(:,end) that differ by less than tol are considered to
%   be equal, and all members of clusters of values where nearest neighbors
%   differ by less than tol are also considered equal.  Note that if tol =
%   0, then no duplicates are eliminated.

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

sordid = sortrows([A; B]);
lastcoldif = abs(diff(sordid(:,end))) >= tol;
othercoldif = any(diff(sordid(:,1:end-1))~=0, 2);
anydif = lastcoldif | othercoldif;
anydif = [true; anydif];  % First element is always a member of unique list.
result = sordid(anydif,:);
