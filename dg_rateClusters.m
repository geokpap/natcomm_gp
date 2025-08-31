function [LRatio, numspikes, clustIDs, LRatioThresh, spikeratio] = ...
    dg_rateClusters(fn1, varargin)
% Reads cluster feature data from .MAT file <fn1>, and returns a list of
% quality ratings with one element for each cluster in the file.  Cluster 0
% is not included in <LRatio>, but it must be supplied in <fn1> in order to
% do the computation.  The quality computation is the L(ratio) described in
% the CLUSTER QUALITY section of Schmitzer-Torbert N & Redish AD (2004) J
% Neurophysiol 91:2259-2272. For explanation of Mahalanobis distance, see
% Fall 2005 web notes for Princeton Computer Science 436 at
% http://www.cs.princeton.edu/courses/archive/fall05/cos436/Duda/PR_Mahal/M_metric.htm. 
%INPUTS
% fn1: a numeric array, or a pathname to a .MAT file containing a numeric
%   array.  If it's the name of a .MAT file, it should contain cluster ID
%   numbers in column 2 and spike feature values in subsequent columns
%   (column 1 is ignored).  If it's a numeric array, it should contain
%   cluster ID numbers in column 1 and spike feature values in subsequent
%   columns; file reading is skipped. Warning: computation time goes as the
%   square of the number of features. The file must contain "unassigned" or
%   "noise" samples as cluster #0 in order to get a legitimate assessment
%   of quality.  The cluster numbers in <fn1> do not need to be
%   consecutive; any missing cluster numbers are simply skipped.
%OUTPUTS
% LRatio: a column vector containing LRatio values for each cluster other
%   than the "noise cluster" number 0.  May contain NaN if 'mahal' silently
%   fails and returns NaN, or if there are too few spikes in the cluster to
%   run 'mahal'.  Indexed by cluster ID numbers, so may also contain NaN
%   for cluster numbers that do not exist.
% numspikes: the number of spikes in each cluster.  Indexed by cluster ID
%   numbers, so may also contain NaN for cluster numbers that do not exist.
% clustIDs: a column vector containing the integer cluster IDs as supplied
%   in <fn1>, in the same order as the values in <LRatio>.  Contains only
%   cluster numbers that exist, and includes cluster "0" if present.
% LRatioThresh: rule of thumb threshold for "good", as judged by eye from
%   waveform overlays and scatter plots (see NOTES).  Indexed by cluster ID
%   numbers, so may also contain NaN for cluster numbers that do not exist.
%   Contains <Inf> when there are too few spikes outside the cluster to run
%   'mahal'.
% spikeratio: the ratio of the number of spikes outside the cluster to the
%   number in the cluster, i.e. <sum(~clustselect) / numspikes(clust_id)>.
%   Indexed by cluster ID numbers, so may also contain NaN for cluster
%   numbers that do not exist.
%OPTIONS
% 'verbose' - reports how many features are in <fn1>.
% 'clusters', clusts - <clusts> is a cell vector, each of whose
%   elements is a vector of indices into <fn1> submitted as a numeric
%   array.  When using 'clusters', the first column of <fn1> is NOT taken
%   to be cluster numbers, but rather is taken to be the first feature.
%   The nth cell in <clusts> is taken to represent cluster <n>.  Cluster
%   0 is considered to be any spikes that do not belong to any of the
%   clusters in <clusts>, and clusters are numbered with consecutive
%   integers starting at 1.
%NOTES
% The worst case (or at least a very bad one) for getting an accurate
% measure of isolation is when the cluster is a hollow shell immediately
% surrounded by another hollow shell of noncluster points.  In this case,
% for a reasonable number of dimensions like 3 or 10, the LRatio numerator
% is roughly
%   0.4 * sum(~clustselect)
% This is where the ratio of the number of points outside the cluster (i.e.
% <sum(~clustselect)>) to the number of points inside <sum(clustselect)>
% becomes important.  When that ratio is 1/8, we get a falsely encouraging
% value of around 0.05 for LRatio for the hollow shell construction in
% spite of the fact the cluster is by construction 100% bogus.  The exact
% value varies slightly with the number of dimensions, but as a rule of
% thumb this is good enough for 3 to 10 dimensions. Therefore we want at
% least as many points outside the cluster as inside in order to get a
% useful result.
%   Empirically, the formula used here:
%       LRatioThresh = 1e-3 * spikeratio + 1e-3
% seems to work pretty well for SE data, with <spikeratio> ranging from 0
% to 2500.
%   Per my analysis of 2-Feb-2022 ("GPClustVol5.odt" pp. 130-154),
% <LRatioThresh> is about 5x lower than it really should be for "general"
% unit analysis, when one wants to select units that are at least
% "semi-good" as judged visually.  Depending on whether one wants to be
% more restrictive or more permissive, one can require satisfaction of both
% of the conditions "LRatio < 0.05" and "LRatio < 5 * LRatioThresh"
% (logical "and"), or one can allow the satisfaction of just one of the two
% conditions (logical "or").

%$Rev: 290 $
%$Date: 2022-04-01 21:48:18 -0400 (Fri, 01 Apr 2022) $
%$Author: dgibson $

clustflag = false;
clusts = {};
verboseflag = false;

argnum = 0;
while true
    argnum = argnum + 1;
    if argnum > length(varargin)
        break
    end
    if ~ischar(varargin{argnum})
        continue
    end
    switch varargin{argnum}
        case 'verbose'
            verboseflag = true;
        case 'clusters'
            argnum = argnum + 1;
            clustflag = true;
            clusts = varargin{argnum};
        otherwise
            error('dg_rateClusters:badoption', ...
                'The option %s is not recognized.', ...
                dg_thing2str(varargin{argnum}));
    end
end

if ischar(fn1)
    % Load the data into arrays s:
    S = load(fn1);
    fields = fieldnames(S);
    s = S.(fields{1});
    clear S;
    s(:,1) = [];
elseif isnumeric(fn1)
    s = fn1;
else
    error('dg_rateClusters:fn1', ...
        '<fn1> must be either char or numeric.');
end

if clustflag
    clustIDs = (1:length(clusts))';
    clustnums = zeros(size(fn1, 1), 1);
    for n = 1:length(clustIDs)
        if any(clustnums(clusts{n}))
            warning('dg_rateClusters:clustIDs', ...
                '<clustIDs> contains overlapping clusters.');
        end
        clustnums(clusts{n}) = n;
    end
    if any(clustnums == 0)
        clustIDs = [0; clustIDs];
    end
    s = [clustnums s];
else
    clustIDs = reshape(unique(s(:,1)), [], 1);
end

if ~ismember(0, clustIDs)
    warning('dg_rateClusters:no0', ...
        'There is no cluster #0');
end
Nfeatures = size(s, 2) - 1;
if Nfeatures < 1
    error('dg_rateClusters:nofeat', ...
        'There are no feature data in %s', fn1 );
end
if verboseflag
    fprintf('Rating clusters based on %d features\n', Nfeatures);
end
L = NaN(max(clustIDs) - 1, 1);
LRatio = NaN(max(clustIDs) - 1, 1);
numspikes = NaN(max(clustIDs) - 1, 1);
LRatioThresh = NaN(max(clustIDs) - 1, 1);
spikeratio = NaN(max(clustIDs) - 1, 1);
for clustidx = 1:length(clustIDs)
    clust_id = clustIDs(clustidx);
    if clust_id == 0
        continue
    end
    % <clustselect> contains all data from the current cluster,
    % clust_id:
    clustselect = (s(:,1) == clust_id);
    numspikes(clust_id, 1) = sum(clustselect);
    spikeratio(clust_id,1) = sum(~clustselect) / numspikes(clust_id);
    LRatioThresh = 7e-4 * spikeratio + 2e-3;
    if sum(clustselect) <= size(s, 2) - 1
        % First column of <s> is IDs, hence the "- 1".
        warning('dg_rateClusters:NaN2', ...
            'Too few spikes in cluster %d to compute LRatio.', ...
            clust_id);
        LRatio(clust_id,1) = NaN;
    else
        % I make the leap of faith that the incredibly compact formulas in
        % the Matlab 'mahal' function really do compute the Mahalanobis
        % distances of Y from mean(X) using the the covariance matrix
        % calculated from X.
        L(clust_id) = sum(1 - chi2cdf(...
            mahal(s(~clustselect, 2:end), s(clustselect, 2:end)), ...
            Nfeatures ));
        LRatio(clust_id,1) = L(clust_id) / numspikes(clust_id);
    end
end

