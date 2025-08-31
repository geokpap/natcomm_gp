function wave = lfp_recipISI(spikechannel, varargin)
%wave = lfp_recipISI(spikechannel) For lfp_createWave; reciprocal ISI
% For use with lfp_createWave.  Computes an instantaneous firing rate
% between each pair of spikes in <spikechannel> as 1/ISI where ISI is the
% Inter-Spike Interval.  No special handling for gaps in recording, so rate
% will appear artifactually low across a gap.  The 1/ISI value extends from
% the sample closest to the first spike and to but not including the sample
% closest to the second spike (except in case of the very last spike, where
% the end sample is included).
%OPTIONS
% 'dur' - the wave represents the duration of the ISI instead of its
%   reciprocal.  Useful for things that aren't really spikes.
% 'verbose' - the usual.

%$Rev: 189 $
%$Date: 2010-12-25 20:04:48 -0500 (Sat, 25 Dec 2010) $
%$Author: dgibson $

lfp_declareGlobals;
durflag = false;
verboseflag = false;
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'dur'
            durflag = true;
        case 'verbose'
            verboseflag = true;
        otherwise
            error('lfp_recipISI:badoption', ...
                ['The option "' dg_thing2str(varargin{argnum}) '" is not recognized.'] );
    end
    argnum = argnum + 1;
end

wave = zeros(1, (length(lfp_TimeStamps) * lfp_SamplesPerFrame));

spike1 = lfp_Spikes{spikechannel}(1); %#ok<*USENS>
sample1 = lfp_time2index(spike1);
numspikes = numel(lfp_Spikes{spikechannel});
for spikeindex = 2:numspikes
    spike2 = lfp_Spikes{spikechannel}(spikeindex);
    isi = spike2 - spike1;
    sample2 = lfp_time2index(spike2);
    if durflag
        wave(sample1:sample2) = isi;
    else
        wave(sample1:sample2) = 1/isi;
    end
    spike1 = spike2;
    sample1 = sample2;
    if verboseflag
        fprintf('Done %.0f/%.0f\n', spikeindex, numspikes);
    end
end