function lfp_plotplace(spikes, xchan, ychan)
%lfp_plotplace(spikes, trial, xchan, ychan)
%   Plots the x and y coordinates of each spike.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

lfp_declareGlobals;
sampleidx = lfp_time2index(spikes);
plot(lfp_Samples{xchan}(sampleidx), lfp_Samples{ychan}(sampleidx), '.');

