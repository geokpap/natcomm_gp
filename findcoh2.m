function findcoh2(sessiondir)
% For Naotaka data.
% Initial screening vs. all LFPs, don't print anything, use a
% post-processor to look through log file for interesting unit-LFP pairs.
% "Interesting" can be defined in terms of the "differential summed
% coherence" or just "significance level of cell counts" (both are logged).
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
% Params for tweaking ====================================================
lfp_XLimAll = [-1 5];
lfp_FreqLim = [20 50];
lfp_AlignmentRef = 2;
% END params for tweaking ================================================

lfp_log(sprintf(...
    '\n\tStarting mycohero2 for %s\n\tlfp_XLimAll=%s, lfp_FreqLim=%s, lfp_AlignmentRef=%s, lfp_CLimAll=%s', ...
    sessiondir, mat2str(lfp_XLimAll), mat2str(lfp_FreqLim), ...
    mat2str(lfp_AlignmentRef), mat2str(lfp_CLimAll) ));

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

% Some trials may not have params 6 - 8, so pad them if necessary:
for k=1:length(lfp_TrialParams); tpsizes(k,:)=size(lfp_TrialParams{k}); end
misfits = find(tpsizes(:,2)~=8);
for k = reshape(misfits,1,[]); lfp_TrialParams{k} = zeros(1,8); end
% Select trial type by condition and number of saccades (e.g. RSQ, 4
% saccades) - <condition>, <numSaccades>, <yellow> = params 6, 7, 8
% respectively.
% Condition can be decoded like below
% 1: RSQ
% 2: FSQ
% 0: LSQ
% 5: Random saccade length
% 3: Trick LSQ
% 4: Reward Delay
% 6: No task with reward
% 7: No task no reward
% 20: saccade without off signal
% 21: Gap single saccade
% 22: No last off signal
% 23: flicker or double reward
% 24: no grid
% 17: Fixed reward schedule rewarded
% 19: Fixed reward schedule not rewarded
% 26: No off reward delay
lfp_selectByRule('HasEvent(lfp_AlignmentRef) && HasParams([6 1; 7 4])');
% Select trial type by number of saccades only (e.g. 4 saccades):
% lfp_selectByRule('HasEvent(lfp_AlignmentRef) && HasParams([7 4])');

enabledtrials = lfp_enabledTrials(find(lfp_SelectedTrials));
lfp_log(sprintf('Selected trials %s', ...
    dg_canonicalSeries(enabledtrials) ));
if isempty(enabledtrials)
    return
end
selectedtrials = lfp_SelectedTrials;

% Select units that have at least mincumrate spikes/s accumulated over all
% trials (i.e. if there are 40 trials and mincumrate = 10, the average rate
% must be at least 0.25 Hz), and that have spikes more uniformly
% distributed than 80% in 1/10 trials.
% Also, Bijan's rule of thumb is to require at least 2 spikes per taper per
% coherence window, so we calculate the mincumrate implied by that.

% Params for tweaking ====================================================
moving_win = [1 .5];
nw = 2;
k = 3;
percentage = 80;
trialfrac = 0.25;
% END params for tweaking ================================================

minspikes = k * 2;
mincumrate = minspikes / moving_win(1);

selectedSpikeChannels = [];
selectedSpikeNames = [];
for channel = 1:length(lfp_Spikes)
    counts = lfp_spikeAnalysis('his', [], channel, [], ...
        'maxskew', percentage, trialfrac);
    % Note that comparisons to NaN are always false:
    if sum(counts) > mincumrate * (lfp_XLimAll(2) - lfp_XLimAll(1))
        selectedSpikeChannels(end+1) = channel;
        selectedSpikeNames = [selectedSpikeNames ' ' lfp_SpikeNames{channel}];
    end
end
lfp_log(sprintf('Selected spike channels %s', selectedSpikeNames));

% Finally do the computation for each selected unit
spikewavechannel = length(lfp_Samples) + 1;
numLFPs = length(lfp_ActiveFilenums);
for channel = reshape(selectedSpikeChannels, 1, [])
    channelname = lfp_SpikeNames{channel};
    cscfilenames = {
        'CSC1.NCS', 'CSC2.NCS', 'CSC3.NCS', 'CSC4.NCS', 'CSC5.NCS', 'CSC6.NCS'
        'CSC7.NCS', 'CSC8.NCS', 'CSC9.NCS', 'CSC10.NCS', 'CSC11.NCS', 'CSC12.NCS'
        'CSC13.NCS', 'CSC14.NCS', 'CSC15.NCS', 'CSC16.NCS', 'CSC17.NCS', 'CSC18.NCS'
        'CSC19.NCS', 'CSC20.NCS', 'CSC21.NCS', 'CSC22.NCS', 'CSC23.NCS', 'CSC24.NCS'
        };
    for round = 1:1
        % lfp_read2('preset', sessiondir, [ {'Events.Nev'}, ...
        %     ]);
        % lfp_SelectedTrials = selectedtrials;
        %lfp_add('preset', sessiondir, ...
            %clusterfiles, 'Naotaka Clusters (*.DWD)', false);
        % spikewavechannel = length(lfp_Samples) + 1;
    if length(lfp_Samples) == spikewavechannel
        lfp_createWave(@lfp_spike2wave, channel, ...
            'name', lfp_SpikeNames{channel}, ...
            'replace', spikewavechannel );
    else
        lfp_createWave(@lfp_spike2wave, channel, ...
            'name', lfp_SpikeNames{channel} );
    end
%        lfp_createWave(@lfp_spike2wave, channel, ...
%            'name', lfp_SpikeNames{channel} );
        for lfpnum = 1:numLFPs
            lfp_log(sprintf('Coherogram: %s vs %s', ...
                lfp_FileNames{lfp_ActiveFilenums(lfpnum)}, ...
                lfp_FileNames{spikewavechannel} ));
            try
                lfp_spec('coh', [], ...
                    [lfp_ActiveFilenums(lfpnum) spikewavechannel], ...
                    moving_win, 'rmdc', 'nw', nw, 'k', k, 'logsig', ...
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
end
close all;

