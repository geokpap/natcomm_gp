function wave = lfp_spike2wave(clusternum, varargin)
%wave = lfp_spike2wave(clusternum); intended for use with lfp_createWave.
% Converts the spike train in spike channel <clusternum> to a CSC wave by
% assigning a 1 to each sample where there is a spike and 0 to all other
% samples.
%OPTIONS
% 'multispike' - actually counts the number of spikes at each sample.  This
%   is useful for low sample rates and situations where spike trains from
%   multiple electrodes have been merged into one multiunit channel.
% 'timestamps' - the first argument is a literal list of spike timestamps
%   rather than the number of a spike channel containing the timestamps.
%   Handy for point process data that you won't use lfp_spikeAnalysis to
%   analyze.

%$Rev: 179 $
%$Date: 2010-12-09 19:49:23 -0500 (Thu, 09 Dec 2010) $
%$Author: dgibson $

lfp_declareGlobals;

multispikeflag = false;
timestampsflag = false;
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'multispike'
            multispikeflag = true;
        case 'timestamps'
            timestampsflag = true;
        otherwise
            error('funcname:badoption', ...
                'The option "%s" is not recognized.', ...
                dg_thing2str(varargin{argnum}) );
    end
    argnum = argnum + 1;
end

if timestampsflag
    samplenums = lfp_time2index(clusternum);
else
    if clusternum > length(lfp_Spikes)
        error('lfp_spike2wave:badclust', ...
            'There is no cluster number %d', clusternum);
    end
    samplenums = lfp_time2index(lfp_Spikes{clusternum});
end
    

wave = zeros(1,length(lfp_TimeStamps) * lfp_SamplesPerFrame);
if multispikeflag
    for idx = 1:length(samplenums)
        wave(samplenums(idx)) = wave(samplenums(idx)) + 1;
    end
else
    wave(samplenums) = 1;
end

