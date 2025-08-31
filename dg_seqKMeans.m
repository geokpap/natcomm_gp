function [idx, C, sumD, D] = dg_seqKMeans(match, indel, X, k, varargin)
%Hacked by Dan Gibson, Oct. 2004 from the Matlab function KMEANS.
%   Additional params <match>, <indel> are used only for edit distance. IDX
%   = DG_SEQKMEANS(match, indel, X, K) partitions the points in the N-by-P
%   data matrix X into K clusters.  This partition minimizes the sum, over
%   all clusters, of the within-cluster sums of point-to-cluster-centroid
%   distances.  X is a cell column vector, each of cell of which contains a
%   sequence of non-negative integers.  KMEANS returns an N-by-1 vector IDX
%   containing the cluster indices of each point.  dg_seqKMeans uses edit
%   distances.
%
%   [IDX, C] = DG_SEQKMEANS(match, indel, X, K) returns the K cluster
%   centroid locations in the K-by-P matrix C.  Values in C represent
%   indices into X.
%
%   [IDX, C, SUMD] = DG_SEQKMEANS(match, indel, X, K) returns the
%   within-cluster sums of point-to-centroid distances in the 1-by-K vector
%   sumD.
%
%   [IDX, C, SUMD, D] = DG_SEQKMEANS(match, indel, X, K) returns distances
%   from each point to every centroid in the N-by-K matrix D.
%
%   [ ... ] = DG_SEQKMEANS(match, indel, ..., 'PARAM1',val1, 'PARAM2',val2,
%   ...) allows you to specify optional parameter name/value pairs to
%   control the iterative algorithm used by KMEANS.  Parameters are:
%
%   'start' - Method used to choose initial cluster centroid positions,
%      sometimes known as "seeds".  Choices are:
%                 {'sample'} - Select K observations from X at random
%                  matrix    - A K-by-1 matrix of starting locations.  In
%                              this case, you can pass in [] for K, and
%                              KMEANS infers K from the first dimension of
%                              the matrix.  You can also supply a 3D array,
%                              implying a value for 'Replicates' from the
%                              array's third dimension. Each element is an
%                              index into the data matrix X.
%
%   'replicates' - Number of times to repeat the clustering, each with a
%      new set of initial centroids [ positive integer | {1}]
%
%   'maxiter' - The maximum number of iterations [ positive integer | {100}]
%
%   'emptyact' - Action to take if a cluster loses all of its member
%      observations.  Choices are:
%               {'error'}    - Treat an empty cluster as an error
%                'drop'      - Remove any clusters that become empty, and
%                              set corresponding values in C and D to NaN.
%                'singleton' - Create a new cluster consisting of the one
%                              observation furthest from its centroid.
%
%   'display' - Display level [ 'off' | {'notify'} | 'final' | 'iter' ]
%
%   Example:
%
%       X = [randn(20,2)+ones(20,2); randn(20,2)-ones(20,2)];
%       [cidx, ctrs] = kmeans(X, 2, 'dist','city', 'rep',5, 'disp','final');
%       plot(X(cidx==1,1),X(cidx==1,2),'r.', ...
%            X(cidx==2,1),X(cidx==2,2),'b.', ctrs(:,1),ctrs(:,2),'kx');
%
%   See also LINKAGE, CLUSTERDATA, SILHOUETTE.

%   KMEANS uses a two-phase iterative algorithm to minimize the sum of
%   point-to-centroid distances, summed over all K clusters.  The first
%   phase uses what the literature often describes as "batch" updates,
%   where each iteration consists of reassigning points to their nearest
%   cluster centroid, all at once, followed by recalculation of cluster
%   centroids. This phase may be thought of as providing a fast but
%   potentially only approximate solution as a starting point for the
%   second phase.  The second phase uses what the literature often
%   describes as "on-line" updates, where points are individually
%   reassigned if doing so will reduce the sum of distances, and cluster
%   centroids are recomputed after each reassignment.  Each iteration
%   during this second phase consists of one pass though all the points.
%   KMEANS can converge to a local optimum, which in this case is a
%   partition of points in which moving any single point to a different
%   cluster increases the total sum of distances.  This problem can only be
%   solved by a clever (or lucky, or exhaustive) choice of starting points.
%
% References:
%
%   [1] Seber, G.A.F., Multivariate Observations, Wiley, New York, 1984.
%   [2] Spath, H. (1985) Cluster Dissection and Analysis: Theory, FORTRAN
%       Programs, Examples, translated by J. Goldschmidt, Halsted Press,
%       New York, 226 pp.

%$Rev: 25 $
%$Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $
%$Author: dgibson $

%
if nargin < 2
    error('At least two input arguments required.');
end

% n points in p dimensional space
[n, p] = size(X);
Xsort = []; Xord = [];

distance = 'edit';
start = 'sample';
reps = [];       
maxit = 100;
emptyact = 'drop';
display = 'notify';

for i=1:2:length(varargin)
	switch varargin{i}
	case 'distance'
		distance = varargin{i+1};
	case 'start'
		start = varargin{i+1};
	case 'replicates'
		reps = varargin{i+1};
	case 'maxiter'
		maxit = varargin{i+1};
	case 'emptyact'
		emptyact = varargin{i+1};
	case 'display'
		display = varargin{i+1};
	end
end

% Allocate distance cache table, initialize diagonal
if ~isequal(class(X), 'cell')
    error('data matrix X must be cell array');
end
global distcache;
distcache = repmat(NaN, n, n);
distcache(find(eye(n))) = 0;

if ischar(start)
    startNames = {'uniform','sample','cluster'};
    i = strmatch(lower(start), startNames);
    if length(i) > 1
        error(sprintf('Ambiguous ''start'' parameter value:  %s.', start));
    elseif isempty(i)
        error(sprintf('Unknown ''start'' parameter value:  %s.', start));
    elseif isempty(k)
        error('You must specify the number of clusters, K.');
    end
    start = startNames{i};
    if strcmp(start, 'uniform')
        if strcmp(distance, 'hamming')
            error('Hamming distance cannot be initialized with uniform random values.');
        end
        Xmins = min(X,1);
        Xmaxs = max(X,1);
    end
elseif isnumeric(start)
    CC = reshape(start, [], 1);
    start = 'numeric';
    if isempty(k)
        k = size(CC,1);
    elseif k ~= size(CC,1);
        error('The ''start'' matrix must have K rows.');
    elseif size(CC,2) ~= p
        error('The ''start'' matrix must have the same number of columns as X.');
    end
    if isempty(reps)
        reps = size(CC,3);
    elseif reps ~= size(CC,3);
        error('The third dimension of the ''start'' array must match the ''replicates'' parameter value.');
    end
    
    % Need to center explicit starting points for 'correlation'. (Re)normalization
    % for 'cosine'/'correlation' is done at each iteration.
    if isequal(distance, 'correlation')
        CC = CC - repmat(mean(CC,2),[1,p,1]);
    end
else
    error('The ''start'' parameter value must be a string or a numeric matrix or array.');
end

if ischar(emptyact)
    emptyactNames = {'error','drop','singleton'};
    i = strmatch(lower(emptyact), emptyactNames);
    if length(i) > 1
        error(sprintf('Ambiguous ''emptyaction'' parameter value:  %s.', emptyact));
    elseif isempty(i)
        error(sprintf('Unknown ''emptyaction'' parameter value:  %s.', emptyact));
    end
    emptyact = emptyactNames{i};
else
    error('The ''emptyaction'' parameter value must be a string.');
end

if ischar(display)
    i = strmatch(lower(display), strvcat('off','notify','final','iter'));
    if length(i) > 1
        error(sprintf('Ambiguous ''display'' parameter value:  %s.', display));
    elseif isempty(i)
        error(sprintf('Unknown ''display'' parameter value:  %s.', display));
    end
    display = i-1;
else
    error('The ''display'' parameter value must be a string.');
end

if k == 1
    error('The number of clusters must be greater than 1.');
elseif n < k
    error('X must have more rows than the number of clusters.');
end

% Assume one replicate
if isempty(reps)
    reps = 1;
end

%
% Done with input argument processing, begin clustering
%

dispfmt = '%6d\t%6d\t%8d\t%12g';
D = repmat(NaN,n,k);   % point-to-cluster distances
Del = repmat(NaN,n,k); % reassignment criterion
m = zeros(k,1);
clustD = repmat(NaN,k,1);    % total distance in cluster for 'edit'

if isequal(distance, 'edit') && ~ismember(start, {'sample', 'numeric'})
    error('That ''start'' option is not available for ''edit'' distance');
end

totsumDBest = Inf;
for rep = 1:reps
    switch start
    case 'sample'
        C = randsample(n,k)'; % contains row indices into X
    case {'numeric'}
        C = CC(:,:,rep);
    end    
    changed = 1:k; % everything is newly assigned
    idx = zeros(n,1);
    totsumD = Inf;
    
    if display > 2 % 'iter'
        disp(sprintf('  iter\t phase\t     num\t         sum'));
    end
    
    %
    % Begin phase one:  batch reassignments
    %
    
    converged = false;
    iter = 0;
    while true
        disp(sprintf(...
            'Phase 1 iter %d (%d max)', iter, maxit ));
        % Compute the distance from every point to each cluster centroid
        D(:,changed) = distfun(X, C(changed), distance, iter, match, indel);
        
        % Compute the total sum of distances for the current configuration.
        % Can't do it first time through, there's no configuration yet.
        if iter > 0
            totsumD = sum(D((idx-1)*n + (1:n)'));
            % Test for a cycle: if objective is not decreased, back out
            % the last step and move on to the single update phase
            if prevtotsumD <= totsumD
                idx = previdx;
                [C(changed), m(changed), clustD(changed)] = gcentroids(X, idx, changed, distance, Xsort, Xord, match, indel);
                iter = iter - 1;
                break;
            end
            if display > 2 % 'iter'
                disp(sprintf(dispfmt,iter,1,length(moved),totsumD));
            end
            if iter >= maxit, break; end
        end

        % Determine closest cluster for each point and reassign points to clusters
        previdx = idx;
        prevtotsumD = totsumD;
        [d, nidx] = min(D, [], 2);

        if iter == 0
            % Every point moved, every cluster will need an update
            moved = 1:n;
            idx = nidx;
            changed = 1:k;
        else
            % Determine which points moved
            moved = find(nidx ~= previdx);
            if length(moved) > 0
                % Resolve ties in favor of not moving
                moved = moved(D((previdx(moved)-1)*n + moved) > d(moved));
            end
            if length(moved) == 0
                break;
            end
            idx(moved) = nidx(moved);

            % Find clusters that gained or lost members
            changed = unique([idx(moved); previdx(moved)])';
        end

        % Calculate the new cluster centroids and counts.
        [C(changed), m(changed), clustD(changed)] = gcentroids(X, idx, changed, distance, Xsort, Xord, match, indel);
        iter = iter + 1;
        
        % Deal with clusters that have just lost all their members
        empties = changed(m(changed) == 0);
        if ~isempty(empties)
            switch emptyact
            case 'error'
                error(sprintf('Empty cluster created at iteration %d.',iter));
            case 'drop'
                % Remove the empty cluster from any further processing
                D(:,empties) = NaN;
                changed = changed(m(changed) > 0);
                if display > 0
                    warning(sprintf('Empty cluster created at iteration %d.',iter));
                end
            case 'singleton'
                if display > 0
                    warning(sprintf('Empty cluster created at iteration %d.',iter));
                end
                
                for i = empties
                    % Find the point furthest away from its current cluster.
                    % Take that point out of its cluster and use it to create
                    % a new singleton cluster to replace the empty one.
                    [dlarge, lonely] = max(d);
                    from = idx(lonely); % taking from this cluster
                    C(i) = lonely;
                    m(i) = 1;
                    idx(lonely) = i;
                    d(lonely) = 0;
                    
                    % Update clusters from which points are taken
                    [C(from), m(from), clustD(from)] = gcentroids(X, idx, from, distance, Xsort, Xord, match, indel);
                    changed = unique([changed from]);
                end
            end
        end
    end % phase one

    % Initialize some information prior to phase two
    
    % index into X specifying the point that is the centroid of each
    % proposed new cluster:
    newCidx = zeros(n,k);   
    newclustD = repmat(NaN,n,k);   
    
    % Ring buffer for detecting limit cycles:
    memlength = min(maxit, 50);
    idxmem = zeros(size(idx,1), memlength);
    
    %
    % Begin phase two:  single reassignments
    %
    changed = find(m' > 0);
    lastmoved = 0;
    nummoved = 0;
    iter1 = iter;
    while iter < maxit
        disp(sprintf(...
            'Phase 2 iter %d (%d max): calc dist', iter, maxit ));
        % Calculate distances to each cluster from each point, and the
        % potential change in total sum of errors for adding or removing
        % each point from each cluster.  Clusters that have not changed
        % membership need not be updated.
        %
        % Singleton clusters are a special case for the sum of dists
        % calculation.  Removing their only point is never best, so the
        % reassignment criterion had better guarantee that a singleton
        % point will stay in its own cluster.  Happily, we get
        % Del(i,idx(i)) == 0 automatically for them.
        previter = iter;
        disp(sprintf('Clusters to update: %s', dg_thing2str(changed)));
        for i = changed
            disp(sprintf('Updating cluster %d', i));
            mbrs = (idx == i);
            mbridx = find(mbrs);
            % For members, compute total distance with the point
            % removed and subtract from original cluster distance:
            for j = 1:length(mbridx);
                newmbridx = mbridx;
                newmbridx(j) = [];
                if isempty(newmbridx)
                    % Theoretically, this move will never actually be
                    % chosen; assigning 0 to newCidx should cause an
                    % error later on if it is.
                    Del(mbridx(j),i) = 0;
                    newCidx(mbridx(j),i) = 0;
                else
                    [newix, newdist] = seqCentroid(X, newmbridx, match, indel);
                    newCidx(mbridx(j),i) = newix;
                    newclustD(mbridx(j),i) = newdist;
                    Del(mbridx(j),i) = clustD(i) - newdist;
                end
            end
            % For nonmembers, add the point and subtract the original
            % cluster distance from the new value
            nonmbridx = find(~mbrs);
            for j = 1:length(nonmbridx);
                newmbridx = mbridx;
                newmbridx(end+1) = nonmbridx(j);
                [newix, newdist] = seqCentroid(X, newmbridx, match, indel);
                newCidx(nonmbridx(j),i) = newix;
                newclustD(nonmbridx(j),i) = newdist;
                Del(nonmbridx(j),i) =  newdist - clustD(i);
            end
        end

        % Determine best possible move, if any, for each point.  Next we
        % will pick one from those that actually did move.
        disp(sprintf(...
            'Phase 2 iter %d (%d max): find moves', iter, maxit ));
        previdx = idx;
        prevtotsumD = totsumD;
        [minDel, nidx] = min(Del, [], 2);
        moved = find(previdx ~= nidx);
        if length(moved) > 0
            % Resolve ties in favor of not moving
            % Apparently, this means to eliminate any moves that don't
            % actually decrease the total distance (-DG)
            moved = moved(Del((previdx(moved)-1)*n + moved) > minDel(moved));
        end
        if length(moved) == 0
            % Count an iteration if phase 2 did nothing at all, or if we're
            % in the middle of FIRST (-DG) pass through all the points
            if (iter - iter1) == 0 | nummoved > 0
                iter = iter + 1;
                if display > 2 % 'iter'
                    disp(sprintf(dispfmt,iter,2,nummoved,totsumD));
                end
            end
            converged = true;
            break;
        end
        
        % Pick the next move in cyclic order
        moved = mod(min(mod(moved - lastmoved - 1, n) + lastmoved), n) + 1;
        
        % If we've gone once through all the points, that's an iteration
        if moved <= lastmoved
            iter = iter + 1;
            if display > 2 % 'iter'
                disp(sprintf(dispfmt,iter,2,nummoved,totsumD));
            end
            if iter >= maxit, break; end
            nummoved = 0;
        end
        nummoved = nummoved + 1;
        lastmoved = moved;
        
        oidx = idx(moved);
        nidx = nidx(moved);
        totsumD = totsumD + Del(moved,nidx) - Del(moved,oidx);
        disp(sprintf('Moved %d from %d to %d', moved, oidx, nidx));
        
        % Update the cluster index vector, and rhe old and new cluster
        % counts and centroids
        idx(moved) = nidx;
        m(nidx) = m(nidx) + 1;
        m(oidx) = m(oidx) - 1;
        switch distance
        case 'edit'
            % Presumably, moved now points at the one point that has moved
            % to a new cluster.
            C(nidx) = newCidx(moved,nidx);
            C(oidx) = newCidx(moved,oidx);
            clustD(nidx) = newclustD(moved,nidx);
            clustD(oidx) = newclustD(moved,oidx);
        end
        changed = sort([oidx nidx]);
        % if end of iteration, check for limit cycle & save current idx
        if iter ~= previter
            foundLimitCycle = false;
            for memiter = 1:memlength
                if all(idx == idxmem(:,memiter))
                    warning(sprintf(...
                        'Found limit cycle on iteration %d of length %d.', ...
                        iter, mod(iter - memiter, memlength) + 1 ));
                    foundLimitCycle = true;
                    break
                end
            end
            if foundLimitCycle == true
                break
            else
                disp(sprintf('iteration: %d', iter));
            end
            idxmem(:, mod(iter, memlength) + 1) = idx;
        end
    end % phase two
    
    if (~converged) & (display > 0)
        warning(sprintf('Failed to converge in %d iterations.', maxit));
    end

    % Calculate cluster-wise sums of distances
    nonempties = find(m(:)'>0);
    D(:,nonempties) = distfun(X, C(nonempties), distance, iter, match, indel);
    d = D((idx-1)*n + (1:n)');
    sumD = zeros(k,1);
    for i = 1:k
        sumD(i) = sum(d(idx == i));
    end
    if display > 1 % 'final' or 'iter'
        disp(sprintf('%d iterations, total sum of distances = %g',iter,totsumD));
    end

    % Save the best solution so far
    if totsumD < totsumDBest
        totsumDBest = totsumD;
        idxBest = idx;
        Cbest = C;
        sumDBest = sumD;
        if nargout > 3
            Dbest = D;
        end
    end
end

% Return the best solution
idx = idxBest;
C = Cbest;
sumD = sumDBest;
if nargout > 3
    D = Dbest;
end


%------------------------------------------------------------------

function D = distfun(X, C, dist, iter, match, indel)
%DISTFUN Calculate point to cluster centroid distances.  There must be a
%column returned for each column in <C>.
% <match> and <indel> are edit penalties for dist = 'edit' and do not have
% to be given otherwise.
[n,p] = size(X);
D = zeros(n,size(C,1));
clusts = 1:size(C,1);

switch dist
case 'edit'
    for i = clusts
        for j = 1:n
            D(j,i) = getdistance(X, j, C(i), match, indel);
        end
    end
end


%------------------------------------------------------------------

function [centroids, counts, sumdist] = gcentroids(X, index, clusts, dist, Xsort, Xord, match, indel)
%GCENTROIDS Centroids and counts (and summed distance for dist = 'edit')
%stratified by group.
[n,p] = size(X);
num = length(clusts);
sumdist = zeros(num, 1);
centroids = zeros(num, 1);
counts = zeros(num,1);
for i = 1:num
    members = find(index == clusts(i));
    if length(members) > 0
        counts(i) = length(members);
        switch dist
        case 'edit'
            seqs = X(members);
            [Cidx, sumdist(i)] = seqCentroid(X, members, match, indel);
            centroids(i) = Cidx;
        end
    end
end

function y = randsample(n, k)
%RANDSAMPLE Random sampling, without replacement
%   Y = RANDSAMPLE(N,K) returns K values sampled at random, without
%   replacement, from the integers 1:N.

%   Copyright 1993-2002 The MathWorks, Inc.
%   $Revision: 25 $  $Date: 2009-03-31 21:56:57 -0400 (Tue, 31 Mar 2009) $

% RANDSAMPLE does not (yet) implement weighted sampling.

if nargin < 2
    error('Requires two input arguments.');
end

% If the sample is a sizeable fraction of the population, just
% randomize the whole population (which involves a full sort
% of n random values), and take the first k.
if 4*k > n
    rp = randperm(n);
    y = rp(1:k);

% If the sample is a small fraction of the population, a full
% sort is wasteful.  Repeatedly sample with replacement until
% there are k unique values.
else
    x = zeros(1,n); % flags
    sumx = 0;
    while sumx < k
        x(ceil(n * rand(1,k-sumx))) = 1; % sample w/replacement
        sumx = sum(x); % count how many unique elements so far
    end
    y = find(x > 0);
    y = y(randperm(k));
end



function [minidx, mindist] = seqCentroid(X, seqs, match, indel)
%SEQCENTROID finds the centermost member of a set of sequences.
%   <seqs> is a vector of row indices into X.  Return the index into X of
%   the sequence with the smallest total distance together with the total
%   distance. (min arbitrarily resolves ties.)

sumdist = zeros(1,length(seqs));

for seqidx = 1:length(seqs)
    others = reshape(seqs,1,[]);
    others(seqidx) = [];
    seq1 = seqs(seqidx);
    for seq2 = others
        sumdist(seqidx) = sumdist(seqidx) + ...
            getdistance(X, seq1, seq2, match, indel);
    end
end

[mindist, minidx] = min(sumdist);
minidx = seqs(minidx);



function d = getdistance(X, seqidx, seqidx2, match, indel)
% <seqidx> and <seqidx2> are row indices into X; calculates value only if it is not
% already stored in distcache.
% This implementation maintains values in lower half-matrix of rectangular
% array distcache.
global distcache;
if seqidx < seqidx2
    i1 = seqidx;
    i2 = seqidx2;
else
    i1 = seqidx2;
    i2 = seqidx;
end
if isnan(distcache(i1,i2))
    distcache(i1,i2) = dg_edDist(X{i1}, X{i2}, match, indel);
end
d = distcache(i1,i2);

    