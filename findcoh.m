function findcoh(sessiondir)
% For Naotaka data.
% Initial screening vs. all LFPs, don't print anything, use a
% post-processor to look through log file for interesting unit-LFP pairs
% defined in terms of the "differential summed coherence", i.e. the sum of
% the magnitude of coherence in all cells in the top half of the
% time-windowed, frequency-windowed, 'signifcolor' coherogram minus that in
% the bottom half, divided by the number of cells in each half. 
%
% Select trials that have lfp_AlignmentRef. Select units subject to minimum
% spike count and maximum skew. Run coherograms against all available LFPs
% for each unit, but do not print anything.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

lfp_declareGlobals;
lfp_changeSetup('naotaka');

% Set the name of the events file here:
lfp_read2('preset', sessiondir, {'Events.Nev'});

% Set time window, frequency window, and alignment event:
lfp_XLimAll = [-1 1];
lfp_FreqLim = [2 20];
lfp_AlignmentRef = 67;

lfp_log(sprintf(...
    '\n\tStarting mycohero2 for %s\n\tlfp_XLimAll=%s, lfp_FreqLim=%s, lfp_AlignmentRef=%s, lfp_CLimAll=%s', ...
    sessiondir, mat2str(lfp_XLimAll), mat2str(lfp_FreqLim), ...
    mat2str(lfp_AlignmentRef), mat2str(lfp_CLimAll) ));
numLFPs = length(lfp_FileNames);

% Find all the *.DWD files for this session:
clusterfiles = {};
files = dir(sessiondir);
for f = files'
    [pathstr,name,ext,versn] = fileparts(upper(f.name));
    if isequal(upper(ext), '.DWD')
        if length(name) < 3
            warning('findcoh:badname1', ...
                'Bad DWD file name: %s', f.name );
        else
            trode = str2num(name(3:end));
            if isempty(trode)
                warning('findcoh:badname2', ...
                    'Bad DWD file name: %s', f.name );
            else
                clusterfiles{end+1} = f.name;
            end
        end
    end
end

if isempty(clusterfiles)
    lfp_log('There are no cluster files in this session');
    return
end
lfp_add('preset', sessiondir, ...
    clusterfiles, 'Naotaka Clusters (*.DWD)', false);

% Select trial type by condition and number of saccades (e.g. RSQ, 4
% saccades):
% lfp_SelectByRule('HasEvent(lfp_AlignmentRef) && HasParams([6 1; 7 4])');
% Select trial type by number of saccades only (e.g. 4 saccades):
lfp_selectByRule('HasEvent(lfp_AlignmentRef) && HasParams([7 4])');

lfp_log(sprintf('Selected trials %s', dg_canonicalSeries(find(lfp_SelectedTrials))));
if ~any(lfp_SelectedTrials)
    return
end

% Select units that have at least mincount spikes/s accumulated over all
% trials (i.e. if there are 40 trials and mincount = 10, the average rate
% must be at least 0.25 Hz), and that have spikes more uniformly
% distributed than 80% in 1/4 trials (this is to get rid of noise bursts,
% units that mysteriously "came and went", etc; the parameters have to be
% adjusted so as NOT to reject units responding to known heterogeneities in
% the task events).  
moving_win = [.5 .1];
mincount = 10;
percentage = 80;
trialfrac = 0.25;
minspikes = mincount * moving_win(1);

filenumpairs = [];
selectedSpikeChannels = [];
selectedSpikeNames = [];
for channel = 1:length(lfp_Spikes)
    counts = lfp_SpikeAnalysis('his', [], channel, [], ...
        'maxskew', percentage, trialfrac);
    % Note that comparisons to NaN are always false:
    if sum(counts) > mincount * (lfp_XLimAll(2) - lfp_XLimAll(1))
        selectedSpikeChannels(end+1) = channel;
        selectedSpikeNames = [selectedSpikeNames ' ' lfp_SpikeNames{channel}];
    end
end
lfp_log(sprintf('Selected spike channels %s', selectedSpikeNames));

% Finally do the computation for each selected unit, re-using the same
% channel for the sampled spike representation each time:
spikewavechannel = length(lfp_Samples) + 1;
for channel = reshape(selectedSpikeChannels, 1, [])
    channelname = lfp_SpikeNames{channel};
    if length(lfp_Samples) == spikewavechannel
        lfp_createWave(@lfp_spike2wave, channel, ...
            'name', lfp_SpikeNames{channel}, ...
            'replace', spikewavechannel );
    else
        lfp_createWave(@lfp_spike2wave, channel, ...
            'name', lfp_SpikeNames{channel} );
    end
    for lfpnum = 1:numLFPs
        lfp_log(sprintf('Coherogram: %s vs %s', ...
            lfp_FileNames{lfp_ActiveFilenums(lfpnum)}, ...
            lfp_FileNames{spikewavechannel} ));
        try
            lfp_spec('coh', [], ...
                [lfp_ActiveFilenums(lfpnum) spikewavechannel], moving_win, ...
                'avg', 'minspikes', minspikes);
        catch
            lfp_log(sprintf(...
                'Error computing coherogram\n%s %s %s %s %s', ...
                dg_thing2str([lfp_ActiveFilenums(lfpnum) spikewavechannel]), ...
                dg_thing2str(moving_win), ...
                'avg', 'minspikes', dg_thing2str(minspikes) ));
        end
    end
end
close all;

