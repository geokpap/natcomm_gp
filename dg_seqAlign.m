function [sim, dist] = dg_seqAlign(v, w, match, indel)
%DEG_SEQALIGN dynamic programming algorithm for sequence alignment
%[sim, dist] = dg_seqAlign(v, w, match, indel)

%[sim, dist] = dg_seqAlign(v, w, match, indel)
% <v> and <w> are arrays of nonnegative integers.  <sim>, <dist> are the
% similarity score and normalized distance respectively between <v> and
% <w>. <match> is the array of incremental similarity scores for matching
% pairs of integers; specifically (because Matlab does not provide
% zero-based indexing), match(p+1,q+1) is the similarity score for matching
% a p with a q. <match> must satisfy match(x,y)=deltamax if x=y,
% match(x,y)<deltamax if x~=y, and match(x,y) = match(y,x) (otherwise <dist>
% may not be a distance measure).  <indel> is a vector of similarities for
% inserting or deleting each integer value; specifically indel(k+1) is the
% incremental similarity score for inserting or deleting the value k.  All
% values in <indel> must be < deltamax.
%
% The values in <indel> are customarily negative, as are the off-diagonal
% elements of <match>.  The diagonal of <match> is customarily 1.  There is
% no special significance to the sign of the scores, however; what really
% matters is the difference between the maximum possible score (deltamax)
% and the score actually assigned to the edit operation.  To compute the
% Longest Common Subsequence, set the off-diagonal elements of <match> =
% -Inf and set indel(:) = 0.
%
% In case of a tie, the (arbitrary) order of preference is del, ins, match.
%
% The 'oops' warning turns out to be a major convoluted issue when
% substitutions are allowed, because the alignment that minimizes distance
% is not necessarily the same as the one that maximizes similarity.

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

if ~isequal(match, match')
    error('dg_seqAlign:asymmetricMatch', ...
        '<match> must be symmetric around the main diagonal.' );
end
deltamax = max(max(match));
if ~isequal(match == deltamax, eye(size(match)))
    error('dg_seqAlign:badDiagMatch', ...
        '<match> must have its maximum value along the main diagonal.' );
end
if (length(indel) ~= size(match,1)) || (length(indel) ~= numel(indel))
    error('dg_seqAlign:badIndel', ...
        '<indel> must be a vector of length equal to one side of <match>.' );
end
if any(indel >= deltamax)
    error('dg_seqAlign:bigIndel', ...
        'All values in <indel> must be < deltamax.' );
end
if ~isequal(fix(v), v) || ~isequal(fix(w), w) || any(v < 0) || any(w < 0)
    error('dg_seqAlign:nonintegers', ...
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
s = zeros(numel(v)+1, numel(w)+1);  % cumulative similarities
d = zeros(numel(v)+1, numel(w)+1);  % cumulative distances
p = zeros(numel(v)+1, numel(w)+1);  % ptrs to predecessors: 1=del, 2=ins, 3=match
p2 = zeros(numel(v)+1, numel(w)+1);  % ptrs to predecessors: 1=del, 2=ins, 3=match
d(1,2:end) = (cumsum(deltamax-indel(w)));
d(2:end,1) = (cumsum(deltamax-indel(v)))';

% Fill in non-boundaries
for i = 2:numel(v)+1
    for j = 2:numel(w)+1
        [s(i,j), p(i,j)] = max(...
            [ s(i-1,j) + indel(v(i-1)), ...
                s(i,j-1) + indel(w(j-1)), ...
                s(i-1,j-1) + match(v(i-1),w(j-1)) ]);
        [d(i,j), p2(i,j)] = min(...
            [ d(i-1,j) + (deltamax - indel(v(i-1))), ...
                d(i,j-1) + (deltamax - indel(w(j-1))), ...
                d(i-1,j-1) + (deltamax - match(v(i-1),w(j-1))) ]);
        if (p2(i,j) ~= p(i,j))
%         if (s(p2) ~= s(p(i,j))) || (d(p2) ~= d(p(i,j)))
            warning('dg_seqAlign:oops', ...
                'Oops - s & d disagree at %d, %d', i, j);
        end
    end
end

sim = s(end,end);
dist = d(end,end);