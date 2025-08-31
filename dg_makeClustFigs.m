function data = dg_makeClustFigs(filename)
%INPUTS
% filename: cut cluster *.mat file with timestamps in col 1, cluster ID
%   numbers in col 2 (0 = "unclustered"), waveforms in cols 3:end.
%OUTPUTS
% data: the array read from <filename>.

%$Rev: 301 $
%$Date: 2022-09-22 18:56:54 -0400 (Thu, 22 Sep 2022) $
%$Author: dgibson $

x = load(filename);
n = fieldnames(x);
data = x.(n{1});
dg_plotClusterWaves(data(:,2), data(:,3:end));
[~, score] = pca(data(:, 3:end));
clustIDs = unique(data(:,2));
clusts = {};
if ismember(0, clustIDs)
    startidx = 2;
else
    startidx = 1;
end
for k = startidx:length(clustIDs)
    clusts{1, end+1} = find(data(:,2) == clustIDs(k)); %#ok<AGROW>
end
hF = dg_plotAllClusterFeatures(clusts, score(:, 1:3), 'size', 4);
hA = subplot(2,2,4,'Parent',hF,'Color','k');
pkval = [max(data(:,3:end), [], 2) min(data(:,3:end), [], 2)];
dg_plotClusterFeatures(clusts, pkval, 'size', 4, 'handle', hA);
xlabel(hA, 'Peak');
ylabel(hA, 'Valley');
dg_plotClusterWavesOverlay(data(:,2), data(:,3:end), ...
    'decimate', round(size(data,1) / 5000));

