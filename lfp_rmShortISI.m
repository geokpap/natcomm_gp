function ts = lfp_rmShortISI(spikechannel, startstop, minISI)
%lfp_getUniformSpikes - for use with lfp_createEvents or 
% lfp_spikeAnalysis(..., 'trigfunc', funcHandle, args).
%
%ts = lfp_getUniformSpikes(spikechannel, startstop, binwidth)
% Samples spikes from spikechannel by removing the second of any pair of spikes
% that are separated by less than <minISI>.  This is done iteratively, so that
% if there are e.g. three spikes in a row s.t. both ISIs are less than
% <minISI>, but the sum of the two ISIs is greater than <minISI>, then only
% the second spike is removed, not the third.
% <minISI> should be set a little higher
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
include = logical(ones(numel(spikes),1));
lastincludedspike = spikes(1);
for k = 2:length(include)
    if spikes(k) - lastincludedspike < minISI
        include(k) = false;
    else
        lastincludedspike = spikes(k);
    end
end
ts = spikes(include);

