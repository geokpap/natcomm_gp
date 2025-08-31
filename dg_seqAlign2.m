function [sim, pairs] = dg_seqAlign2(v, w, match, indel)
%DEG_SEQALIGN2 dynamic programming algorithm for sequence alignment
%[sim, pairs] = dg_seqAlign2(v, w, match, indel)

%[sim, pairs] = dg_seqAlign2(v, w, match, indel)
% Modifed from dg_seqAlign2 by DG 1/2/2007.
% <v> and <w> are arrays of nonnegative integers.  <sim> is the accumulated
% similarity score between <v> and <w>.  <pairs> is a two-column array
% containing the indices of matched elements in <v> and <w>, respectively.
% <match> is the array of incremental similarity scores for matching [sim,
% pairs] of integers; specifically (because Matlab does not provide
% zero-based indexing), match(p+1,q+1) is the similarity score for matching
% a p with a q. <match> must satisfy match(x,y)=deltamax if x=y,
% match(x,y)<deltamax if x~=y, and match(x,y) = match(y,x) (otherwise
% <dist> may not be a distance measure).  <indel> is a vector of
% similarities for inserting or deleting each integer value; specifically
% indel(k+1) is the incremental similarity score for inserting or deleting
% the value k.  All values in <indel> must be < deltamax.
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

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

if numel(v) * max(max(match)) > 65535 || numel(w) * max(max(match)) > 65535
    % Addition to a uint16 silently clips at 65535:
    error('dg_seqAlign2:overflow', ...
        '<v> and <w> must each be fewer than (65536/max(max(match))) elements.' );
end
if numel(v) * numel(w) > 0.7e+9
    % Practically speaking, on a win32 box, even the nojvm version will run
    % out of memory (2 GB for Matlab + data combined) somewhere around
    % here:
    error('dg_seqAlign2:overflow2', ...
        'numel(v) * numel(w) must be < 0.7e+9.' );
end
if ~isequal(match, match')
    error('dg_seqAlign2:asymmetricMatch', ...
        '<match> must be symmetric around the main diagonal.' );
end
deltamax = max(max(match));
if ~isequal(match == deltamax, eye(size(match)))
    error('dg_seqAlign2:badDiagMatch', ...
        '<match> must have its maximum value along the main diagonal.' );
end
if (length(indel) ~= size(match,1)) || (length(indel) ~= numel(indel))
    error('dg_seqAlign2:badIndel', ...
        '<indel> must be a vector of length equal to one side of <match>.' );
end
if any(indel >= deltamax)
    error('dg_seqAlign2:bigIndel', ...
        'All values in <indel> must be < deltamax.' );
end
if ~isequal(fix(v), v) || ~isequal(fix(w), w) || any(v < 0) || any(w < 0)
    error('dg_seqAlign2:nonintegers', ...
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
s = uint16(zeros(numel(v)+1, numel(w)+1));  % cumulative similarities
p = int8(zeros(numel(v)+1, numel(w)+1));  % ptrs to predecessors: 1=del, 2=ins, 3=match

% Fill in non-boundaries
for i = 2:numel(v)+1
    for j = 2:numel(w)+1
        [s(i,j), p(i,j)] = max(...
            [ s(i-1,j) + indel(v(i-1)), ...
                s(i,j-1) + indel(w(j-1)), ...
                s(i-1,j-1) + match(v(i-1),w(j-1)) ]);
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
        error('dg_seqAlign2:badp', ...
            'Internal error.' );
    end
end
if pairnum
    pairs(1:pairnum, :) = [];
end

