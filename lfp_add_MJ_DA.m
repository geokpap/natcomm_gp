function lfp_add_MJ_DA(dirname)
% lfp_add_MJ_DA(dirname)
%   Unpacks DA data from MJ-format.mat files into CSC channels.  A new CSC
%   channel is created for each channel named by the fields in
%   meanDAData.allTrial.
%INPUT
% dirname: optional arg to specify directory from which to load files.
%   If not given, the current value of lfp_DataDir is used.
%   The data files must be named rawdata.mat and dadata.mat, and must
%   contain the following variables.
%   rawdata.mat:
%       Evt: struct with field 'EvtTS' that contains Neuralynx
%           timestamps of (some? all?) events.
%       EvtPos: struct with field 'TTLout' that contains indices into
%           Evt.EvtTS for the Nlx/voltammetry sync events.  In particular,
%           EvtPos.trialStart(trial) is the index of the trial start event.
%           Synchronization between the Nlx clock and the voltammetry clock
%           is accomplished by means of a synchronization TTL event s.t.
%           for each trial, Evt.EvtTS(EvtPos.TTLout) matches 
%           TTLs.(posTTL,1).
%       TTLs: 3-column array where col. 1 contains voltammetry system
%           timestamps for each scan, in seconds.
%       posTTL: struct with field 'Bit241' that contains indices into
%           TTLs for the Nlx/voltammetry sync events.
%   dadata.mat:
%       meanDAData: struct with field 'allTrial' which is a struct with
%           one field for each voltammetry channel, each of which contains
%           a cell array with one element per trial where each element is
%           a row vector of DA coefficients for every scan in the trial.
%           The names of the fields in meanDAData.allTrial are used as
%           channel names in lfp_FileNames.

global lfp_DataDir lfp_SamplePeriod lfp_Samples lfp_TrialIndex
scanperiod = 0.1; % seconds between scans

if nargin < 1
    dirname = lfp_DataDir;
end

rawdata = load(fullfile(dirname, 'rawdata.mat'));
dadata = load(fullfile(dirname, 'dadata.mat'));
chnames = fieldnames(dadata.meanDAData.allTrial);
numch = length(chnames);
numtrials = length(dadata.meanDAData.allTrial.(chnames{1}));
if numtrials ~= size(lfp_TrialIndex,1)
    error('lfp_add_MJ_DA:trials', ...
        'The number of trials in %s does not match the number in lfp_TrialIndex', ...
        fullfile(dirname, 'dadata.mat'));
end
for chidx = 2:length(chnames)
    if length(dadata.meanDAData.allTrial.(chnames{chidx})) ~= numtrials
        error('lfp_add_MJ_DA:meanDAData', ...
            'meanDAData contains different numbers of trials for different channels.');
    end
end
for chidx = 1:numch
    if isempty(lfp_Samples)
        error('lfp_add_MJ_DA:samples', ...
            'There must already be at least one CSC channel in memory.');
    end
    interpdata = NaN(size(lfp_Samples{1}));
    lastsamp = 0;
    for trial = 1:size(lfp_TrialIndex,1)
        trialdata = dadata.meanDAData.allTrial.(chnames{chidx}){trial};
        trialstartTS = rawdata.Evt.EvtTS(rawdata.EvtPos.trialStart(trial));
        scan1TS = trialstartTS - 1; % -10 * 0.1 = -1 s re: trial start
        scan1samp = lfp_time2index(scan1TS);
        trialscanTS = (0 : length(trialdata) - 1) * scanperiod + ...
            scan1TS;
        numtrialsamp = ceil( (length(trialdata) - 1) * scanperiod / ...
            lfp_SamplePeriod );
        trialsampTS = lfp_index2time(scan1samp) + ...
            (0 : (numtrialsamp - 1)) * lfp_SamplePeriod;
        if lastsamp ~= 0 && lastsamp < scan1samp - 1
            % Interpolate between trials:
            interpdata(lastsamp : scan1samp - 1) = interp1( ...
                [lastsamp scan1samp], ...
                [interpdata(lastsamp) trialdata(1)], ...
                lastsamp : scan1samp - 1 );
            if scan1samp - lastsamp - 1 > ceil( ...
                    scanperiod / lfp_SamplePeriod )
                warning('lfp_add_MJ_DA:gap', ...
                    'There is a gap of %d s between trials %d and %d', ...
                    (scan1samp - lastsamp) * lfp_SamplePeriod, ...
                    trial - 1, trial );
            end
        elseif chidx == 1 && lastsamp ~= 0 && lastsamp > scan1samp
            warning('lfp_add_MJ_DA:overlap', ...
                    'There is an overlap of %d s between trials %d and %d', ...
                    (lastsamp - scan1samp) * lfp_SamplePeriod, ...
                    trial - 1, trial );
        end
        lastsamp = scan1samp + numtrialsamp - 1;
        % Interpolate between scans of current trial:
        interpdata(scan1samp + (0 : (numtrialsamp - 1))) = interp1( ...
            trialscanTS, trialdata, trialsampTS, 'linear', 'extrap' );
    end
    lfp_createWave(@lfp_waverecord, 1, interpdata, 'name', ...
        chnames{chidx}, 'units', 'nM');
end
