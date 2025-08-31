function [bitvolts, ADmaxval, numdisorderedTS, maxes, mins, dcs, meds, ...
    fractionClipped] = dg_histCSCStats(CSCStats)
% Creates a figure with histograms of the values in CSCStats aggregated
% across all files.
%INPUT
% CSCStats: as saved by dg_computeCSCStats, i.e. as loaded by
%	load(fullfile(sessiondir, 'computeCSCStats.mat'));
%OUTPUTS
% bitvolts, ADmaxval, numdisorderedTS, maxes, mins, dcs, meds: values
%   aggregated across all files in <CSCStats>, e.g. <meds> contains all
%   values of CSCStats.medianpeak.

bitvolts = cell2mat({CSCStats.bitvolts});
ADmaxval = cell2mat({CSCStats.ADmaxval});
numdisorderedTS = cellfun(@length, {CSCStats.disorderedTSidx});
maxes = cell2mat({CSCStats.max});
mins = cell2mat({CSCStats.min});
dcs = cell2mat({CSCStats.dc});
meds = cell2mat({CSCStats.medianpeak});
fractionClipped = cell2mat({CSCStats.fractionClipped});

figure;
subplot(4, 2, 1);
if any(~isnan(bitvolts))
    hist(bitvolts, linspace(0, max(bitvolts), 100));
end
title('bitvolts');
subplot(4, 2, 2);
if any(~isnan(ADmaxval))
    hist(ADmaxval, linspace(0, max(ADmaxval), 100));
end
title('ADmaxval');
subplot(4, 2, 3);
hist(numdisorderedTS, linspace(0, max(numdisorderedTS)+1, 100));
title('numdisorderedTS');
subplot(4, 2, 4);
hist(maxes);
title('maxes');
subplot(4, 2, 5);
hist(mins);
title('mins');
subplot(4, 2, 6);
hist(dcs);
title('dcs');
subplot(4, 2, 7);
hist(meds);
title('median peaks');
subplot(4, 2, 8);
hist(fractionClipped);
title('fractionClipped');

