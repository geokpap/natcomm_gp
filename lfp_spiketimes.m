function ts = lfp_spiketimes(clusternum)
% LFP_SPIKETIMES is just a dumb wrapper to make spike times available to
% lfp_createEvents.
%ts = lfp_spiketimes(clusternum)

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

lfp_declareGlobals;

if clusternum > length(lfp_Spikes)
    error('lfp_spiketimes:badclust', ...
        'There is no cluster number %d', clusternum);
end
ts = lfp_Spikes{clusternum};