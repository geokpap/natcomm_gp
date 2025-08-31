function isosc = dg_IsoSc(clustdata, lambda)
%Compute the "Isolation Score" of Joshua, Elias, Levine, Bergman; Journal of
%Neuroscience Methods 163 (2007) 267?282.
%INPUTS
% clustdata: cluster ID numbers in column 1 and feature values in
%   subsequent columns.  One row per observation.
% lambda: gain for Similarity measure.
%OUTPUTS
% isosc: a column vector containing Iso.Sc. values for each cluster other
%   than the "noise cluster" number 0.  Indexed by cluster ID numbers, so
%   may contain NaN for cluster numbers that do not exist.

%$Rev:  $
%$Date:  $
%$Author: dgibson $

clustID = reshape(unique(clustdata(:,1)), [], 1);
isosc = NaN(max(clustID), 1);
for clustidx = 1:length(clustID)
    if clustID(clustidx) == 0
        continue
    end
    isinclust = clustdata(:,1) == clustID(clustidx);
    % Per eqn 8, X and Y are two individual points, and d(X, Y) is
    % the Euclidean distance between them.  <d0> is the average Euclidian
    % distance in the spike cluster <clustidx>.  <lambda> is a free
    % parameter to control the gain of the Similarity =
    % exp(-d(X,Y)*lambda/d0).  <PxY> is a scalar value for the two
    % individual points in eqn 9.  To make the notation readable, we
    % represent the the sets X, Y, and Z as vectors of row indices into
    % <clustdata>, which we designate simply as
    % <X>, <Y>, and <Z>.
    X = reshape(find(isinclust), [], 1);
    numpts = length(X);
    
    %
    % Compute <d0>.
    %
    dists = NaN(numpts);
    Y = X;
    % We avoid duplicating the same distance calculation by computing only
    % the upper half-matrix, and copying to lower half.
    for xidx = 1:numpts-1
        dists(xidx, xidx) = 0;
        for yidx = xidx+1:numpts
            dists(xidx, yidx) = d(X(xidx), Y(yidx), clustdata);
            dists(yidx, xidx) = dists(xidx, yidx);
        end
    end
    % It's unclear whether they meant to include the self-distances of 0 or
    % not, but for consistency with the treatment of <Z> (below) we exclude
    % them here.
    d0 = mean(dists(~eye(numpts)));
    
    %
    % Compute the normalized similarity matrix <PxY>.
    %
    PxY = NaN(numpts);
    for xidx = 1:numpts
        for yidx = 1:numpts
            Z = [1:(X(xidx)-1) (X(xidx)+1):size(clustdata,1)];
            denom = sum( ...
                exp(-d(repmat(X(xidx), 1, length(Z)), Z, clustdata) ...
                * lambda / d0 ) );
            if yidx ~= xidx
                PxY(xidx, yidx) = ...
                    (exp(-dists(xidx, yidx) * lambda / d0)) / denom;
            end
        end
    end
    % The text introducing eqn 8 states, "we compute the normalized
    % similarity between each event in the spike cluster, X, to all other
    % events (spikes and noise), Y". Taken literally, that implies that the
    % self-similarities of 1 are to be excluded.  However, since the
    % self-term is explicitly excluded from the denominator in eqn (9), I
    % take the liberty of excluding it here as well, as far as the sums are
    % concerned, by setting the diagonal terms to 0.
    PxY(logical(eye(numpts))) = 0;
    if isequal(PxY, PxY')
        warning('dg_IsoSc:PxY', 'PxY shows commutativity again.');
    end

    %
    % Calculate <isosc>.
    %
    % Here, <Px> is a vector of values for all points <X> in the cluster.
    Px = sum(PxY,2); % eqn 10; sum is over all <Y> in cluster.
    % Per eqn 11, this sum is over all points <X> in cluster:
    isosc(clustID(clustidx),1) = sum(Px) / numpts; % eqn 11
end

end


function result = d(X, Y, clustdata)
% Euclidean distance between points <X> and <Y>.  <X> and <Y> must
% be of the same size.
%INPUTS
% X: numerical row index into <clustdata>.
% Y: numerical row index into <clustdata>.
%OUTPUT
% result: array of distances from each element of <X> to the
%   corresponding element of <Y>; same size as <X> and <Y>.
if ~isequal(size(X), size(Y))
    error('dg_IsoSc:d', ...
        '<X> and <Y> must be same size.');
end
vecdiff = clustdata(X,2:3) - clustdata(Y,2:3);
result = reshape( sqrt(sum(vecdiff.^2, 2)), size(X) );
end
