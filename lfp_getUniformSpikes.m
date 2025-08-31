function ts = lfp_getUniformSpikes(spikechannel, startstop, binwidth)
%lfp_getUniformSpikes - for use with lfp_createEvents or 
% lfp_spikeAnalysis(..., 'trigfunc', funcHandle, args).
%
%ts = lfp_getUniformSpikes(spikechannel, startstop, binwidth)
% Samples spikes from spikechannel by taking the first spike in each disjoint
% time bin <binwidth> seconds wide.  <binwidth> should be set a little higher
% than the "typical" (mean, median, whatever) ISI to eliminate 
% over-representation of times of high firing rate, but note that times of
% low firing rate will still be under-represented in that case.  <startstop>
% is the two-element vector of indices into lfp_Events required by 'trigfunc';
% if it is empty, then all of the time covered by the spike channel is used.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

lfp_declareGlobals;

spikes = lfp_Spikes{spikechannel};
if ~isempty(startstop)
    spikes = spikes( ...
        spikes >= lfp_Events(startstop(1),1)...
        & spikes < lfp_Events(startstop(2),1) );
end
numbins = ceil((spikes(end) - spikes(1))/binwidth);
ts = NaN(numbins,1);
binedges = spikes(1) + (0:numbins)*binwidth;
for k = 1:numbins
    inbin = find( spikes >= binedges(k)...
        & spikes < binedges(k+1) );
    if ~isempty(inbin)
        ts(k) = spikes(inbin(1));
    end
end
ts(isnan(ts)) = [];

