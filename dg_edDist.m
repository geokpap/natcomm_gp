function [dist, ptr] = dg_edDist(v, w, match, indel)
%DG_EDDIST dynamic programming algorithm for sequence alignment
%[dist, ptr] = dg_edDist(v, w, match, indel)

%[dist, ptr] = dg_edDist(v, w, match, indel)
% <v> and <w> are arrays of nonnegative integers.  <dist> is the  distance
% between <v> and <w>. 
% <ptr> is the pointer matrix used to calculate <dist>, where 1=del (up),
% 2=ins (left), 3=match (diagonal).
%
% <match> is the array of incremental distance scores for matching pairs of
% integers; specifically (because Matlab does not provide zero-based
% indexing), match(p+1,q+1) is the similarity score for matching a p with a
% q. <match> must satisfy match(x,y)=0 if x=y, match(x,y)>0 if x~=y, and
% match(x,y) = match(y,x) (otherwise <dist> may not be a distance measure).
% <indel> is a vector of distances for inserting or deleting each integer
% value; specifically indel(k+1) is the incremental distance score for
% inserting or deleting the value k.  All values in <indel> must be >0.
%
% To compute the Longest Common Subsequence, set the off-diagonal elements
% of <match> = Inf and set indel(:) = 1.
%
% In case of a tie, the (arbitrary) order of preference is del, ins, match.

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

if ~isequal(match, match')
    error('dg_edDist:asymmetricMatch', ...
        '<match> must be symmetric around the main diagonal.' );
end
if ~isequal(match == 0, eye(size(match)))
    error('dg_edDist:badDiagMatch', ...
        '<match> must be zero along the main diagonal and only there.' );
end
if (length(indel) ~= size(match,1)) || (length(indel) ~= numel(indel))
    error('dg_edDist:badIndel', ...
        '<indel> must be a vector of length equal to one side of <match>.' );
end
if any(indel < 0)
    error('dg_edDist:bigIndel', ...
        'All values in <indel> must be >0.' );
end
if ~isequal(fix(v), v) || ~isequal(fix(w), w) || any(v < 0) || any(w < 0)
    error('dg_edDist:nonintegers', ...
        '<w> and <v> must be nonnegative integer arrays.' );
end
v = reshape(v,[],1);
w = reshape(w,1,[]);
indel = reshape(indel, 1, []);

%Matlabize the sequences (i.e. make them strictly positive integers that
%can be used to index match and indel directly)
w = w + 1;
v = v + 1;

% Initialize arrays
d = zeros(numel(v)+1, numel(w)+1);  % cumulative distances
ptr = zeros(numel(v)+1, numel(w)+1);  % ptrs to predecessors: 1=del, 2=ins, 3=match
d(1,2:end) = (cumsum(indel(w)));
d(2:end,1) = (cumsum(indel(v)))';

% Fill in non-boundaries
for i = 2:numel(v)+1
    for j = 2:numel(w)+1
        i1 = i - 1;
        j1 = j - 1;
        [d(i,j), ptr(i,j)] = min(...
            [ d(i1,j) + indel(v(i1)), ...
                d(i,j1) + indel(w(j1)), ...
                d(i1,j1) + match(v(i1),w(j1)) ]);
    end
end

dist = d(end,end);