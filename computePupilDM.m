function computePupilDM(sessionID)
% Computes and saves decision matrix plots of average pupil diameter during
% the entire cue period, separately for ApAv and ApAp.

global lfp_TrialParams lfp_AlignmentRef lfp_Samples

% Find <sessiondir> from <sessionID>
dataroot = '/annex2/analysis/dgibson/eyelink3';
cmdstr = sprintf('find %s -type d -name %s', dataroot, sessionID);
[status, result] = system(cmdstr);
if status
    error('computePupilDM:nosession', ...
        'Failed to find session directory.');
end
resultlines = strsplit(result, char(10)); %#ok<CHARTEN>
resultlines(cellfun(@isempty, resultlines)) = [];
if length(resultlines) ~= 1
    error('computePupilDM:resultlines', ...
        '%d resultlines', length(resultlines));
end
sessiondir = resultlines{1};
fprintf('sessiondir: %s\n', sessiondir);

% Find and read required files
lfp_changeSetup('georgios_ApAvApAp');
evtnames = {'Events.mat' 'events.mat'};
pupilnames = {'csc126clean_down32.mat'};
evtfile = '';
for k = 1:length(evtnames)
    if exist(fullfile(sessiondir, evtnames{k}), 'file')
        evtfile = evtnames{k};
        break
    end
end
if isempty(evtfile)
    error('computePupilDM:noEvt', ...
        'No events file.');
end
pupilfile = '';
for k = 1:length(pupilnames)
    if exist(fullfile(sessiondir, pupilnames{k}), 'file')
        pupilfile = pupilnames{k};
        break
    end
end
if isempty(evtfile)
    error('computePupilDM:noEvt', ...
        'No pupil file.');
end
lfp_read2('preset', sessiondir, {evtfile, pupilfile});
lfp_getEvtIDs;
blinklev = -0.05;
lfp_Samples{1}(lfp_Samples{1}(:) < blinklev) = NaN;

% Select completed choice trials for ApAv:
lfp_AlignmentRef = choiceCueOn; %#ok<*NASGU>
lfp_selectByRule( ...
    'HasEvent(lfp_AlignmentRef) && HasEvent([airpuffOnAp rewardOnAv])');
[trials, filenums, window, ~, ...
    ~, getSamplesOpts] = lfp_CSCboilerplate({[], 1, [0 1.5]});
[sampledata, ~, ~, ~, badtrials, trials, ...
    ~, ~] = lfp_getSamples( trials, filenums, window, ...
    getSamplesOpts{:} );
trialsApAv = setdiff(trials, badtrials);
paramsApAv = cell2mat(lfp_TrialParams(trialsApAv, :));
avgdiamApAv = nanmean(sampledata, 1);
[resultApAv, binedgesX, binedgesY] = dg_smooth2Dscatter( ...
    paramsApAv(:, [2 1]), avgdiamApAv', 'bounds', [0, 190, 0, 190], ...
    'numbins', 19 ); %#ok<ASGLU>

% Select completed choice trials for ApAp:
lfp_AlignmentRef = choiceCueOnApAp;
lfp_selectByRule( ...
    'HasEvent(lfp_AlignmentRef) && HasEvent([rewardOnRedApAp rewardOnYelApAp])');
[trials, filenums, window, ~, ...
    ~, getSamplesOpts] = lfp_CSCboilerplate({[], 1, [0 1.5]});
[sampledata, ~, ~, ~, badtrials, trials, ...
    ~, ~] = lfp_getSamples( trials, filenums, window, ...
    getSamplesOpts{:} );
trialsApAv = setdiff(trials, badtrials);
paramsApAv = cell2mat(lfp_TrialParams(trialsApAv, :));
avgdiamApAv = nanmean(sampledata, 1);
[resultApAp, binedgesX, binedgesY] = dg_smooth2Dscatter( ...
    paramsApAv(:, [2 1]), avgdiamApAv', 'bounds', [0, 190, 0, 190], ...
    'numbins', 19 );

% Wrap-up
resultdiff = resultApAv - resultApAp;
save( fullfile(sessiondir, 'computePupilDM.mat'), 'resultdiff', ...
    'resultApAp', 'resultApAv', 'binedgesX', 'binedgesY', 'sessiondir' );
fprintf('Done %s\n', sessiondir);
