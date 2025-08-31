function [TS, clustIDs, Samples, bitvolts, retrig, scores, PCs, frac] = ...
    dg_mergeNoiseCluster(SEfile, clusterfile)
% Read raw Neuralynx spike snapshot file and file containing cluster
% assignments for each time stamp assigned to a cluster.  Match
% assigned time stamps to raw time stamps, extract waveform exclusion
% criteria ("invalidation" criteria) from assigned clusters, remove all
% waveforms that are invalidated by those criteria, and return spike data
% for assigned and valid unassigned waveforms only, with cluster
% assignments, including cluster 0 for unassigned waveforms.
%INPUTS
% SEfile: absolute or relative pathname for Neuralynx spike snapshot file
%   to use.
% clusterfile: absolute or relative pathname for a text file containing an
%   ASCII rendition of the 'Multiple Cut Cluster (*.MAT)' format.
%   File should contain array with spike timestamps in seconds in col
%   1 and cluster ID in col 2 as first variable saved.
%OUTPUTS
% TS: timestamps in microseconds for all "valid" waveforms.
% clustIDs: integer cluster IDs, size(TS).  Includes 0 to designate valid
%   waveforms that were not assigned to any cluster.
% Samples: as returned by dg_readSpike.
% bitvolts: from <SEfile> header; multiply <Samples> by <bitvolts> to get
%   sample values in volts.
% scores: projections of each spike's waveform on each principal component,
%   with one spike per row, and the projections arranged in columns.
% PCs: the principal components on which each spike gets projected, with
%   one principal component per column.
% frac: fraction of total merged spikes that were invalidated.
%NOTES
% Currently only works for *.nse files as <SEfile>. 

%$Rev: 274 $
%$Date: 2021-08-02 14:17:09 -0400 (Mon, 02 Aug 2021) $
%$Author: dgibson $

if exist(SEfile, 'file')
    [TS, Samples, header] = dg_readSpike(SEfile);
else
    error('dg_mergeNoiseCluster:SEfile', ...
        'SEfile "%s" does not exist.', SEfile);
end
for k = 1:length(header)
    if regexp(header{k}, '^\s*-ADBitVolts\s+')
        ADBitVoltstr = regexprep(header{k}, ...
            '^\s*-ADBitVolts\s+', '');
        bitvolts = str2double(ADBitVoltstr);
        break
    else
        bitvolts = NaN;
    end
end
for k = 1:length(header)
    if regexp(header{k}, '^\s*-SpikeRetriggerTime\s+')
        retrigstr = regexprep(header{k}, ...
            '^\s*-SpikeRetriggerTime\s+', '');
        retrig = str2double(retrigstr);
        break
    else
        retrig = NaN;
    end
end
clustIDs = zeros(size(TS)); % cluster numbers for each element of <TS>
clusts = load(clusterfile, '-ascii');
if isempty(clusts)
    error('dg_mergeNoiseCluster:empty', ...
        'Clusterfile %s is empty.', clusterfile);
end
% All we need is timestamps and cluster numbers:
if size(clusts,2) > 2
    clusts(:,3:end) = [];
end
% Convert <clusts> to microseconds:
clusts(:,1) = round(clusts(:,1) * 1e6);

% Iterate through the assigned <clusts> and mark assignments in <clustIDs>:
TSidx = 1;
for wfidx = 1:size(clusts, 1)
    while TSidx <= length(TS) && clusts(wfidx, 1) > TS(TSidx)
        TSidx = TSidx + 1;
    end
    if TS(TSidx) ~= clusts(wfidx, 1)
        error('dg_mergeNoiseCluster:missing', ...
            'The spike at %d in clusterfile does not exist in SEfile', ...
            clusts(wfidx, 1));
    end
    clustIDs(TSidx) = clusts(wfidx, 2);
end

% Find upper and lower envelopes of cluster-assigned waveforms, as row
% vectors, and use them to define invalidation lines.  Invalidation lines
% are calculated as the envelope of the top and bottom 0.1% of sample
% values from waveforms, expanded by a factor of 1.1 from the mean
% envelope.
startpt = 9; % starting sample for invalidations
clustenvs = prctile( squeeze(Samples(startpt:end, 1, clustIDs ~= 0))', ...
    [0.1 99.9] );
meanenv = mean(clustenvs);
envdiff = clustenvs(2,:) - clustenvs(1,:);
invalidations(1,:) = meanenv + 1.1 * envdiff / 2;
invalidations(2,:) = meanenv - 1.1 * envdiff / 2;
isinvalid = false(size(TS));
for wfidx = find(clustIDs == 0)
    isinvalid(wfidx) = ...
        any(Samples(startpt:end, 1, wfidx)' > invalidations(1,:)) ...
        || any(Samples(startpt:end, 1, wfidx)' < invalidations(2,:));
end
if all(isinvalid(clustIDs == 0))
    error('dg_mergeNoiseCluster:allinvalid', ...
        'All waveforms in cluster 0 were invalidated.');
end

frac = sum(isinvalid) / length(TS);
TS(isinvalid) = [];
clustIDs(isinvalid) = [];
Samples(:, :, isinvalid) = [];

[PCs, scores] = pca(squeeze(Samples)');

 