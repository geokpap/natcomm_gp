function spikes = lfp_getTrialSpikes(trial, clustnums, preRoll, postRoll)
%   Returns a cell array containing a column vector of spike times that
%   belong to <trial> for each spike channel in <clustnums>.  The time
%   interval belonging to <trial> is extended in the negative direction by
%   <preRoll> and in the positive direction by <postRoll> when deciding
%   which spikes to include.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

lfp_declareGlobals;

timerange = [ lfp_Events(lfp_TrialIndex(trial,1), 1) - preRoll ...
    lfp_Events(lfp_TrialIndex(trial,2), 1) + postRoll ];
for clustidx = 1:length(clustnums)
    spikes{clustidx}(:,1) = lfp_Spikes{clustnums(clustidx)}(find( ...
        lfp_Spikes{clustnums(clustidx)} > timerange(1) ...
        & lfp_Spikes{clustnums(clustidx)} < timerange(2) )) ;
end