function [spikewaves, spikecounts, clustIDs] = dg_lookupSpikeWaveSE( ...
    clusts, TS, samp)
% Functionally analogous to dg_lookupSpikeWaveforms, but for Single
% Electrode data.  Also, this does not attempt to do any kind of fancy
% stuff with finding or even reading files, it just operates on
% already-read data.  If there is a cluster #0 in <clusts>, it is ignored.
% Any spikes in <TS> that do not belong to a cluster with a nonzero ID
% number are considered to belong to cluster #0
%INPUTS
% clusts: a numeric array containing timestamps in the first column,
%   cluster ID numbers in the second column, and one row for each spike.
%   May contain additional columns, but they are ignored.
% TS: Timestamps read from raw SE (*.nse) file.
% samp: Spike waveform samples read from raw SE (*.nse) file; one row per
%   spike, one column per sample (note this is not the format returned by
%   dg_readSpike).
%OUTPUTS
% spikewaves: one row per spike; col. 1 = timestamp, col. 2 = clustIDs, col.
%   3:end = waveform samples.
% spikecounts: one row per value of clustIDs; col. 1 = number of spikes in
%   <clusts>, col. 2 = number of spikes in <TS> and therefore in
%   <spikewaves>.
% clustIDs: column vector of numerical cluster IDs as given in <clusts>, in
%   ascending order, not including cluster #0 if there is one.
%OPTIONS
%NOTES

%$Rev:  $
%$Date:  $
%$Author: dgibson $

tol = 1.01e-4;
clustIDs = reshape(unique(clusts(:,2)), [], 1);
clustIDs(clustIDs == 0) = [];
TS = reshape(TS, [], 1);
spikewaves = [TS zeros(length(TS), 1) samp];
spikecounts = NaN(length(clustIDs), 2);
for clustidx = 1:length(clustIDs)
    fprintf('dg_lookupSpikeWaveSE clustidx=%d\n', clustidx);
    isInClust = dg_tolIsmember( TS, ...
        clusts(clusts(:,2) == clustIDs(clustidx), 1), tol );
    spikecounts(clustidx, 1) = sum(clusts(:,2) == clustIDs(clustidx));
    spikecounts(clustidx, 2) = sum(isInClust);
    if spikecounts(clustidx, 1) ~= spikecounts(clustidx, 2)
        warning('dg_lookupSpikeWaveSE:spikecount', ...
            'For clustIDs=%d, <clusts> contains %d and <TS> contains %d spikes.', ...
            clustIDs(clustidx), sum(clusts(:,2) == clustIDs(clustidx)), ...
            sum(isInClust) );
    end
    spikewaves(isInClust, 2) = clustIDs(clustidx);
end

