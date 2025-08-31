function bulkProcess_bill(processnum, dirlist, varargin)
%BULKPROCESS_BILL performs offline lfp_lib processing on a list of
%directories.
% bulkProcess_bill(processnum, dirlist)
%   <processnum> selects from a collection of predefined processing
%   protocols.  <dirlist> is a cell array containing the list of session
%   directories that are to be processed; if only one directory is to be
%   processed, it must still be a cell array, not a string.

%   PROCESSING PROTOCOLS
%       1: mycoherograms - initial screening, printsig
%       2: mydebunker - log shuffles and print coh vs. all 6 LFP only for
%       units listed in lookup table
%       3: scrambledLFPcoherograms - 
%       4: mycohero2 - initial screening vs. all LFPs, don't print
%       anything, use a post-processor to look through log file for
%       interesting unit-LFP pairs
%       5: rmEPcoh - same as mycohero2, except that it removes the EP from
%       both the LFP and the spike wave before computing coherence, and it
%       only processes units listed in lookup table.  NOTE: since minspikes
%       does not work with EP-removed spike waveforms, it is not used here,
%       so be sure that you have enough spikes before running this
%       procedure.
%       6: BPvsEvt - "Band Power vs Events"; for each trode, compute
%       lfp_BandPower for lfp_FreqLim, then  compute the average over
%       trials and over time in a window around each event.  Save to a
%       spreadsheet, 'BPvsEvtLog.xls', the average over trodes of power vs.
%       event, the one trode with the greatest range, and each trode
%       individually.  The actual event IDs are specified in a hard-coded
%       list 'evtIDs', and each has its own value for lfp_XLimAll.
%       Columns: sessiondir, lfp_FreqLim, maxtrode, avg1...N, max1...N.
%       The last two items are the power for each event averaged over
%       electrodes and the power on the one electrode with the greatest
%       range between min and max power as a function of alignment event.
%       Power is raw, NOT normalized.
%       7: makeParthaPack - uses lfp_export to create a MAT-file containing
%       all event, LFP, and wave-converted spike data for a rat session.
%       8: BPvsVelo - construct and save scatter plots of band-limited
%       power (averaged over all electrodes) vs. velocity and acceleration.
%       Saves summary stats (r and p for velocity and accleration) to
%       spreadsheet "BPvsVelo.xls".
%       9: mycohero3 - initial screening of all LFPs vs. all LFPs, don't
%       print anything, use a post-processor to look through log file for
%       interesting unit-LFP pairs
%       10: makeParthaPack2 - same as makeParthaPack, except it sets
%       lfp_SelectedTrials and lfp_XLimAll to avoid gate artifact, and then
%       applies lfp_SelectedTrials and lfp_XlimAll before saving. Settings
%       are chosen as follows:  selects trials that have the
%       lfp_AlignmentRef; then runs lfp_selectByDuration using lfp_XLimAll
%       = [-1 5] for lfp_AlignmentRef = 23 or 14, and lfp_XLimAll = [-4 2]
%       for lfp_AlignmentRef = [15 16] or [17 18]; then looks at each
%       trial's wave data to find the last saturated (2047 or -2048) data
%       point within the time period common to lfp_XLimAll and all the
%       selected trials, then searches trial by trial to find the 3rd
%       latest last saturated point, and rounds in quarter-second
%       increments forwards from that point to set a new left-hand limit of
%       lfp_XLimAll.  Any trials that still contain saturated data points
%       within the common time determined by the new lfp_XLimAll (at most
%       two) are then added to lfp_BadTrials.  Output filenames are of the
%       form <animalID><sessionID>wSp_e<evtID>_v6.mat, where evt_ID is the
%       value of the first element of lfp_AlignmentRef.

% OPTIONS
%   'debug' - simply executes the function without try...catch error
%   handling.

% As a coding convenience, each process is defined by a function that
% accepts a session directory as its only argument.  The functions can be
% named anything, but the name must be explicitly listed as a handle in
% processHandles.  If the order of the function handles is changed in
% processHandles, then the corresponding value of processnum is also
% changed.

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

warning off MATLAB:warn_r14_function_handle_transition
processHandles = [ @mycoherograms
    @mydebunker
    @scrambledLFPcoherograms
    @mycohero2
    @rmEPcoh
    @BPvsEvt
    @makeParthaPack
    @BPvsVelo
    @mycohero3
    @makeParthaPack2
];

debugflag = false;
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'debug'
            debugflag = true;
        otherwise
            error('bulkProcess_bill:badoption', ...
                ['The option "' varargin{argnum} ...
                    '" is not recognized.'] );
    end
    argnum = argnum + 1;
end

if processnum < 1 || processnum > length(processHandles)
    error('bulkProcess_joey:badprocessnum', ...
        '<processnum> must be between 1 and %d inclusive', ...
        length(processHandles) );
end
for diridx = 1:length(dirlist)
    mydir = dirlist{diridx};
    if debugflag
        feval(processHandles(processnum), mydir);
    else
        try
            feval(processHandles(processnum), mydir);
        catch
            lfp_declareGlobals;
            if isempty(lfp_LogFileName)
                logname = 'lfp_lib.log';
                lfp_LogFileName = which(logname);
                if isempty(lfp_LogFileName)
                    lfp_LogFileName = fullfile(pwd, logname);
                end
            end
            [msgstr, msgid] = lasterr;
            logmsg = sprintf('Error while processing %s\n%s\n%s', ...
                mydir, msgid, msgstr );
            lfp_log(logmsg);
            disp(logmsg);
        end
    end
end
warning on MATLAB:warn_r14_function_handle_transition


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function definitions for processing protocols %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mycoherograms(sessiondir)
% Sets lfp_SelectedTrials.
% Using the existing values of all other globals:
% Select trials that have lfp_AlignmentRef.
% Select trials by duration to maximize total analyzed time.
% Select units subject to minimum spike count and maximum skew.
% Print avg coherogram for each selected unit and pop avg for session.
% For each coherogram printed, run 3 seeded shuffle coherograms but do not
% print them (significantish ones will be flagged in lfp_lib.log).
% We assume here that the LFP files in sessiondir are numbered sequentially
% starting at 1, with no skips.

lfp_declareGlobals;

if isempty(lfp_XLimAll) || isempty(lfp_AlignmentRef)
    error('lfp_XLimAll and lfp_AlignmentRef must have non-empty values for ''mycoherograms''');
end

lfp_read2('preset', sessiondir, {'Events.Nev'});
lfp_log(sprintf(...
    '\n\tStarting mycoherograms for %s\n\tlfp_XLimAll=%s, lfp_FreqLim=%s, lfp_AlignmentRef=%s, lfp_CLimAll=%s', ...
    sessiondir, mat2str(lfp_XLimAll), mat2str(lfp_FreqLim), ...
    mat2str(lfp_AlignmentRef), mat2str(lfp_CLimAll) ));
numLFPs = length(lfp_FileNames);

% Find all the *.Tn and *.Tnn files for this session:
clusterfiles = {};
[animaldir, sessionname] = fileparts(sessiondir);
files = dir(animaldir);
for f = files'
    [pathstr,name,ext,versn] = fileparts(upper(f.name));
    if (length(ext) >= 2) && (ext(2) == 'T') ...
            && isequal(upper(sessionname), name(1:length(sessionname)))
        trode = str2num(ext(3:end));
        if ~isempty(trode)
            clusterfiles{end+1} = f.name;
        end
    end
end

lfp_add('preset', animaldir, ...
    clusterfiles, 'Rodent Clusters (*.Tnn, *.TTn)', false);

lfp_selectByRule('HasEvent(lfp_AlignmentRef)');
lfp_selectByDuration('and');
lfp_log(sprintf('Selected trials %s', dg_canonicalSeries(find(lfp_SelectedTrials))));
if ~any(lfp_SelectedTrials)
    return
end

% Select units that have at least mincumrate spikes/s accumulated over all
% trials (i.e. if there are 40 trials and mincumrate = 10, the average rate
% must be at least 0.25 Hz), and that have spikes more uniformly
% distributed than 80% in 1/10 trials:
moving_win = [.5 .1];
mincumrate = 10;
percentage = 80;
trialfrac = 0.1;
minspikes = mincumrate * moving_win(1);

filenumpairs = [];
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
for channel = reshape(selectedSpikeChannels, 1, [])
    channelname = lfp_SpikeNames{channel};
    Tidx = find(channelname == 'T');
    Cidx = find(channelname == 'C');
    trodenum = str2num(channelname(Tidx(end)+1:Cidx(end)-1));
    if trodenum > numLFPs
        warning('bulkProcess_bill:badtrodenum', ...
            'Electrode %d has no LFP', trodenum );
        lfp_log(sprintf('Electrode %d has no LFP', trodenum));
    else
        lfp_createWave(@lfp_spike2wave, channel, 'name', lfp_SpikeNames{channel});
        filenumpairs = [ filenumpairs trodenum lfp_ActiveFilenums(end) ];
        lfp_log(sprintf('Coherogram: %s vs %s', ...
            lfp_FileNames{filenumpairs(end-1)}, ...
            lfp_FileNames{filenumpairs(end)} ));
        lfp_spec('coh', [], filenumpairs(end-1:end), moving_win, 'avg', ...
            'rmdc', 'nw', 2, 'k', 3, ...
            'printsig', 'minspikes', minspikes, 'signifcolor', 'antigray');
%         for seed = 1:3 
%             rand('state', seed);
%             lfp_log(['Set random seed = ', num2str(seed)]);
%             lfp_log(sprintf('Shuffled Coherogram: %s vs %s', ...
%                 lfp_FileNames{filenumpairs(end-1)}, ...
%                 lfp_FileNames{filenumpairs(end)} ));
%             lfp_spec('coh', [], filenumpairs(end-1:end), moving_win, ...
%             'rmdc', 'nw', 2, 'k', 3, ...
%                 'avg', 'shuffle', 'printsig', 'minspikes', minspikes, 'antigray');
%         end
    end
end
% lfp_log(sprintf('Population Coherogram: %s', ...
%     mat2str(filenumpairs) ));
% lfp_spec('coh', [], filenumpairs, moving_win, 'pop', 'printsig', ...
%             'rmdc', 'nw', 2, 'k', 3, ...
%     'minspikes', minspikes, 'signifcolor', 'antigray');
% for seed = 1:3 
%     rand('state', seed);
%     lfp_log(['Set random seed = ', num2str(seed)]);
%     lfp_log(sprintf('Shuffled Pop Coherogram: %s', ...
%         mat2str(filenumpairs) ));
%     lfp_spec('coh', [], filenumpairs, moving_win, 'pop', 'shuffle',
%     'printsig', 'minspikes', minspikes, 'antigray', ...
%         'rmdc', 'nw', 2, 'k', 3, ...
%         );
% end
close all;



function mydebunker(sessiondir)
% Runs the shuffles, does not print (just logs)
% Also prints all 6 non-shuffled coherograms of the unit vs. all 6 LFPs

lfp_declareGlobals;
moving_win = [.5 .1];
mincumrate = 10;
percentage = 80;
trialfrac = 0.1;
minspikes = mincumrate * moving_win(1);

if isempty(lfp_XLimAll) || isempty(lfp_AlignmentRef)
    error('lfp_XLimAll and lfp_AlignmentRef must have non-empty values');
end

lfp_read2('preset', sessiondir, {'Events.Nev'});
lfp_log(sprintf(...
    '\n\tStarting mydebunker for %s\n\tlfp_XLimAll=%s, lfp_FreqLim=%s, lfp_AlignmentRef=%s, lfp_CLimAll=%s', ...
    sessiondir, mat2str(lfp_XLimAll), mat2str(lfp_FreqLim), ...
    mat2str(lfp_AlignmentRef), mat2str(lfp_CLimAll) ));
numLFPs = length(lfp_FileNames);

% Find all the *.Tn and *.Tnn files for this session:
clusterfiles = {};
[animaldir, sessionname] = fileparts(sessiondir);
files = dir(animaldir);
for f = files'
    [pathstr,name,ext,versn] = fileparts(upper(f.name));
    if (length(ext) >= 2) && (ext(2) == 'T') ...
            && isequal(upper(sessionname), name(1:length(sessionname)))
        trode = str2num(ext(3:end));
        if ~isempty(trode)
            clusterfiles{end+1} = f.name;
        end
    end
end

lfp_add('preset', animaldir, ...
    clusterfiles, 'Rodent Clusters (*.Tnn, *.TTn)', false);

lfp_selectByRule('HasEvent(lfp_AlignmentRef)');
lfp_selectByDuration('and');
lfp_log(sprintf('Selected trials %s', dg_canonicalSeries(find(lfp_SelectedTrials))));
if ~any(lfp_SelectedTrials)
    return
end

% Find units from lookup table:
theunits = {
    '\\Matrisome2\RData\Bill_LFP_Analysis\B04\acq05' {'ACQ05-T1C4'}
};

tablerow = find(ismember(theunits(:,1), sessiondir));
if isempty(tablerow)
    error('mydebunker:nosuchdir', ...
        'The session directory %s is not in the lookup table', sessiondir );
end

filenumpairs = [];
selectedSpikeChannels = [];
selectedSpikeNames = [];
for channel = 1:length(lfp_Spikes)
    if ismember(lfp_SpikeNames(channel), theunits{tablerow, 2})
        selectedSpikeChannels(end+1) = channel;
        selectedSpikeNames = [selectedSpikeNames ' ' lfp_SpikeNames{channel}];
    end
end
lfp_log(sprintf('Selected spike channels %s', selectedSpikeNames));
for channel = reshape(selectedSpikeChannels, 1, [])
    channelname = lfp_SpikeNames{channel};
    Tidx = find(channelname == 'T');
    Cidx = find(channelname == 'C');
    trodenum = str2num(channelname(Tidx(end)+1:Cidx(end)-1));
    lfp_createWave(@lfp_spike2wave, channel, 'name', lfp_SpikeNames{channel});
    for seed = 1:3 
        rand('state', seed);
        lfp_log(['Set random seed = ', num2str(seed)]);
        % disambiguate lfp_ActiveFilenums for stupid Matlab interpreter:
        lfp_ActiveFilenums = lfp_ActiveFilenums; 
        lfp_log(sprintf('Shuffled Coherogram: %s vs %s', ...
            lfp_FileNames{trodenum}, ...
            lfp_FileNames{lfp_ActiveFilenums(end)} ));
        lfp_spec('coh', [], [trodenum lfp_ActiveFilenums(end)], moving_win, ...
            'rmdc', 'nw', 2, 'k', 3, ...
            'avg', 'shuffle', 'minspikes', minspikes, 'antigray');
    end
    for trodenum = 1:6
        lfp_log(sprintf('Coherogram: %s vs %s', ...
            lfp_FileNames{trodenum}, ...
            lfp_FileNames{lfp_ActiveFilenums(end)} ));
        lfp_spec('coh', [], [trodenum lfp_ActiveFilenums(end)], moving_win, 'avg', ...
            'rmdc', 'nw', 2, 'k', 3, ...
            'print', 'minspikes', minspikes, 'signifcolor', 'antigray');
    end
end
close all;




function scrambledLFPcoherograms(sessiondir)
% Same as mycoherograms, but uses a lookup table to derive the
% corresponding LFP number from the electrode number.

lfp_declareGlobals;

if isempty(lfp_XLimAll) || isempty(lfp_AlignmentRef)
    error('lfp_XLimAll and lfp_AlignmentRef must have non-empty values');
end

lfp_read2('preset', sessiondir, {'Events.Nev'});
lfp_log(sprintf(...
    '\n\tStarting scrambledLFPcoherograms for %s\n\tlfp_XLimAll=%s, lfp_FreqLim=%s, lfp_AlignmentRef=%s, lfp_CLimAll=%s', ...
    sessiondir, mat2str(lfp_XLimAll), mat2str(lfp_FreqLim), ...
    mat2str(lfp_AlignmentRef), mat2str(lfp_CLimAll) ));
numLFPs = length(lfp_FileNames);

% Find all the *.Tn and *.Tnn files for this session:
clusterfiles = {};
[animaldir, sessionname] = fileparts(sessiondir);
files = dir(animaldir);
for f = files'
    [pathstr,name,ext,versn] = fileparts(upper(f.name));
    if (length(ext) >= 2) && (ext(2) == 'T') ...
            && isequal(upper(sessionname), name(1:length(sessionname)))
        trode = str2num(ext(3:end));
        if ~isempty(trode)
            clusterfiles{end+1} = f.name;
        end
    end
end

lfp_add('preset', animaldir, ...
    clusterfiles, 'Rodent Clusters (*.Tnn, *.TTn)', false);

lfp_selectByRule('HasEvent(lfp_AlignmentRef)');
lfp_selectByDuration('and');
lfp_log(sprintf('Selected trials %s', dg_canonicalSeries(find(lfp_SelectedTrials))));
if ~any(lfp_SelectedTrials)
    return
end

% Select units that have at least mincumrate spikes/s accumulated over all
% trials (i.e. if there are 40 trials and mincumrate = 10, the average rate
% must be at least 0.25 Hz), and that have spikes more uniformly
% distributed than 80% in 1/10 trials:
moving_win = [.5 .1];
mincumrate = 10;
percentage = 80;
trialfrac = 0.1;
minspikes = mincumrate * moving_win(1);

filenumpairs = [];
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
% trodemap maps from electrode number to LFP number (i.e. the number in
% parentheses is the electrode, the one to the right of the '=' is the LFP);
% any electrode that is not listed here will be skipped with a warning.
trodemap(3) = 2;
trodemap(6) = 5;
trodemap(8) = 3;
trodemap(9) = 3;
trodemap(10) = 6;
for channel = reshape(selectedSpikeChannels, 1, [])
    channelname = lfp_SpikeNames{channel};
    Tidx = find(channelname == 'T');
    Cidx = find(channelname == 'C');
    trodenum = str2num(channelname(Tidx(end)+1:Cidx(end)-1));
    if trodenum < length(trodemap)
        lfpnum = trodemap(trodenum);
    else
        lfpnum = 0;
    end
    if lfpnum == 0
        warning('bulkProcess_bill:badlfpnum', ...
            'Electrode %d has no LFP in trodemap', trodenum );
        lfp_log(sprintf('Electrode %d has no LFP', trodenum));
    else
        lfp_createWave(@lfp_spike2wave, channel, 'name', lfp_SpikeNames{channel});
        filenumpairs = [ filenumpairs lfpnum lfp_ActiveFilenums(end) ];
        lfp_log(sprintf('Coherogram: %s vs %s', ...
            lfp_FileNames{filenumpairs(end-1)}, ...
            lfp_FileNames{filenumpairs(end)} ));
        lfp_spec('coh', [], filenumpairs(end-1:end), moving_win, 'avg', ...
            'rmdc', 'nw', 2, 'k', 3, ...
            'printsig', 'minspikes', minspikes, 'signifcolor', 'antigray');
    end
end
close all;



function mycohero4(sessiondir)
% Sets lfp_SelectedTrials.
% Using the existing values of all other globals:
% Select trials that have lfp_AlignmentRef.
% Select trials by duration to maximize total analyzed time.
% Select units subject to minimum spike count and maximum skew.
% Run coherograms against all available LFPs for each unit, but do not
% print anything.
% We assume here that the LFP files in sessiondir are numbered sequentially
% starting at 1, with no skips.

lfp_declareGlobals;

if isempty(lfp_XLimAll) || isempty(lfp_AlignmentRef)
    error('lfp_XLimAll and lfp_AlignmentRef must have non-empty values');
end

lfp_read2('preset', sessiondir, {'Events.Nev'});
lfp_log(sprintf(...
    '\n\tStarting mycohero2 for %s\n\tlfp_XLimAll=%s, lfp_FreqLim=%s, lfp_AlignmentRef=%s, lfp_CLimAll=%s', ...
    sessiondir, mat2str(lfp_XLimAll), mat2str(lfp_FreqLim), ...
    mat2str(lfp_AlignmentRef), mat2str(lfp_CLimAll) ));
numLFPs = length(lfp_ActiveFilenums);

% Find all the *.Tn and *.Tnn files for this session:
clusterfiles = {};
[animaldir, sessionname] = fileparts(sessiondir);
files = dir(animaldir);
for f = files'
    [pathstr,name,ext,versn] = fileparts(upper(f.name));
    if (length(ext) >= 2) && (ext(2) == 'T') ...
            && isequal(upper(sessionname), name(1:length(sessionname)))
        trode = str2num(ext(3:end));
        if ~isempty(trode)
            clusterfiles{end+1} = f.name;
        end
    end
end

lfp_add('preset', animaldir, ...
    clusterfiles, 'Rodent Clusters (*.Tnn, *.TTn)', false);

lfp_selectByRule('HasEvent(lfp_AlignmentRef)');
lfp_selectByDuration('and');
lfp_log(sprintf('Selected trials %s', dg_canonicalSeries(find(lfp_SelectedTrials))));
if ~any(lfp_SelectedTrials)
    return
end

% Select units that have at least mincumrate spikes/s accumulated over all
% trials (i.e. if there are 40 trials and mincumrate = 10, the average rate
% must be at least 0.25 Hz), and that have spikes more uniformly
% distributed than 80% in 1/10 trials.
% Also, Bijan's rule of thumb is to require at least 2 spikes per taper per
% coherence window, so we calculate the mincumrate implied by that.

% Params for tweaking ====================================================
moving_win = [1 .25];
nw = 1.5;
k = 3;
percentage = 80;
trialfrac = 0.1;
% END params for tweaking ================================================

minspikes = k * 2;
mincumrate = minspikes / moving_win(1);

filenumpairs = [];
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
for channel = reshape(selectedSpikeChannels, 1, [])
    channelname = lfp_SpikeNames{channel};
    lfp_createWave(@lfp_spike2wave, channel, 'name', lfp_SpikeNames{channel});
    for lfpnum = 1:numLFPs
        % disambiguate lfp_ActiveFilenums for stupid Matlab interpreter:
        lfp_ActiveFilenums = lfp_ActiveFilenums; 
        lfp_log(sprintf('Coherogram: %s vs %s', ...
            lfp_FileNames{lfpnum}, ...
            lfp_FileNames{lfp_ActiveFilenums(end)} ));
        hF = lfp_spec('coh', [], [lfp_ActiveFilenums(lfpnum) lfp_ActiveFilenums(end)], moving_win, ...
            'rmdc', 'nw', nw, 'k', k, 'logsig', ...
            'avg', 'minspikes', minspikes, 'lin');
        saveas(hF, sprintf('%s%s_%s%s_%d.fig', ...
            animalID, sessionname, ...
            lfp_FileNames{lfpA}, lfp_FileNames{lfpB}, ...
            lfp_AlignmentRef(1) ));
    end
end
close all;




function rmEPcoh(sessiondir)
% Sets lfp_SelectedTrials.
% Using the existing values of all other globals:
% Select trials that have lfp_AlignmentRef.
% Select trials by duration to maximize total analyzed time.
% Select units subject to minimum spike count and maximum skew.
% Run coherograms against all available LFPs for each unit, but do not
% print anything.
% We assume here that the LFP files in sessiondir are numbered sequentially
% starting at 1, with no skips.

lfp_declareGlobals;

if isempty(lfp_XLimAll) || isempty(lfp_AlignmentRef)
    error('lfp_XLimAll and lfp_AlignmentRef must have non-empty values');
end

lfp_read2('preset', sessiondir, {'Events.Nev'});
lfp_log(sprintf(...
    '\n\tStarting rmEPcoh for %s\n\tlfp_XLimAll=%s, lfp_FreqLim=%s, lfp_AlignmentRef=%s, lfp_CLimAll=%s', ...
    sessiondir, mat2str(lfp_XLimAll), mat2str(lfp_FreqLim), ...
    mat2str(lfp_AlignmentRef), mat2str(lfp_CLimAll) ));
numLFPs = length(lfp_ActiveFilenums);

% Find all the *.Tn and *.Tnn files for this session:
clusterfiles = {};
[animaldir, sessionname] = fileparts(sessiondir);
files = dir(animaldir);
for f = files'
    [pathstr,name,ext,versn] = fileparts(upper(f.name));
    if (length(ext) >= 2) && (ext(2) == 'T') ...
            && isequal(upper(sessionname), name(1:length(sessionname)))
        trode = str2num(ext(3:end));
        if ~isempty(trode)
            clusterfiles{end+1} = f.name;
        end
    end
end

lfp_add('preset', animaldir, ...
    clusterfiles, 'Rodent Clusters (*.Tnn, *.TTn)', false);

lfp_selectByRule('HasEvent(lfp_AlignmentRef)');
lfp_selectByDuration('and');
lfp_log(sprintf('Selected trials %s', dg_canonicalSeries(find(lfp_SelectedTrials))));
if ~any(lfp_SelectedTrials)
    return
end

% Select units that have at least mincumrate spikes/s accumulated over all
% trials (i.e. if there are 40 trials and mincumrate = 10, the average rate
% must be at least 0.25 Hz), and that have spikes more uniformly
% distributed than 80% in 1/10 trials.
% Also, Bijan's rule of thumb is to require at least 2 spikes per taper per
% coherence window, so we calculate the mincumrate implied by that.

% Params for tweaking ====================================================
moving_win = [1 .1];
nw = 3;
k = 5;
percentage = 80;
trialfrac = 0.1;
% END params for tweaking ================================================

minspikes = k * 2;
mincumrate = minspikes / moving_win(1);

filenumpairs = [];
selectedSpikeChannels = [];
selectedSpikeNames = [];
% Find units from lookup table:
theunits = {
    'C:\B03\Acq05' {'ACQ05-T6C2'}
    'C:\B08\Acq24' {'ACQ24-T6C3'}
    };

tablerow = find(ismember(theunits(:,1), sessiondir));
if isempty(tablerow)
    error('rmEPcoh:nosuchdir', ...
        'The session directory %s is not in the lookup table', sessiondir );
end

filenumpairs = [];
selectedSpikeChannels = [];
selectedSpikeNames = [];
for channel = 1:length(lfp_Spikes)
    if ismember(lfp_SpikeNames(channel), theunits{tablerow, 2})
        selectedSpikeChannels(end+1) = channel;
        selectedSpikeNames = [selectedSpikeNames ' ' lfp_SpikeNames{channel}];
    end
end
lfp_log(sprintf('Selected spike channels %s', selectedSpikeNames));
for channel = reshape(selectedSpikeChannels, 1, [])
    channelname = lfp_SpikeNames{channel};
    lfp_createWave(@lfp_spike2wave, channel, 'name', lfp_SpikeNames{channel});
    % disambiguate lfp_ActiveFilenums for stupid Matlab interpreter:
    lfp_ActiveFilenums = lfp_ActiveFilenums;
    wavechannel = lfp_ActiveFilenums(end);
    % temporarily enlarge lfp_XLimAll for EP reomval:
    oldxlimall = lfp_XLimAll;
    widerxlim = lfp_XLimAll + [-(moving_win(1)/2+lfp_SamplePeriod) moving_win(1)/2+lfp_SamplePeriod];
    lfp_XLimAll = widerxlim;
    % remove EP from spike wave:
    [hF, data] = lfp_disp([], ...
        [wavechannel], ...
        [], 'avg');
    close(hF);
    warning('off', 'lfp_rmEP:noRef');
    lfp_createWave(@lfp_rmEP, wavechannel, ...
        data{1}(1,1), data{1}(:,2), ...
        'name', [lfp_FileNames{wavechannel} 'rmEP']);
    warning('on', 'lfp_rmEP:noRef');
    newwavechannel = length(lfp_FileNames);
    for lfpnum = 1:numLFPs
        % remove EP from LFP:
        lfp_XLimAll = widerxlim;
        [hF, data] = lfp_disp([], ...
            [lfp_ActiveFilenums(lfpnum)], ...
            [], 'avg');
        close(hF);
        warning('off', 'lfp_rmEP:noRef');
        lfp_createWave(@lfp_rmEP, lfp_ActiveFilenums(lfpnum), ...
            data{1}(1,1), data{1}(:,2), ...
            'name', [lfp_FileNames{lfp_ActiveFilenums(lfpnum)} 'rmEP']);
        warning('on', 'lfp_rmEP:noRef');
        newlfpchannel = length(lfp_FileNames);
        % log a message and compute coherence
        lfp_XLimAll = oldxlimall;
        lfp_log(sprintf('Coherogram: %s vs %s', ...
            lfp_FileNames{newlfpchannel}, ...
            lfp_FileNames{newwavechannel} ));
        lfp_spec('coh', [], ...
            [newlfpchannel newwavechannel], moving_win, ...
            'rmdc', 'nw', nw, 'k', k, 'logsig', ...
            'avg');
    end
end
close all;



function BPvsEvt(sessiondir)
lfp_declareGlobals;
if isempty(lfp_FreqLim)
    error('lfp_FreqLim must have non-empty value');
end
% Specify event IDs, etc %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note that this function is NOT smart enough to know when you change
% these, so if you do change them, you should rename any existing output
% file so that the new results for the new events will go in a new file
% with new header rows.
evtIDs = {[31 38] [31 38] 23 14 [15 16] [17 18]};
xlimalls = { [-1 0] [0 1] [-.5 .5] [0 1] [-.5 .5] [-.5 .5] };
moving_win = [1 .25];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if numel(evtIDs) ~= numel(xlimalls)
    error('bulkProcess_bill:badlen', ...
    'evtIDs and xlimalls are different lengths' );
end
lfp_read2('preset', sessiondir, {'Events.Nev'});
outfilename = fullfile(pwd, 'BPvsEvtLog.xls');
if exist(outfilename) == 2
    newfile = false;
else
    newfile = true;
end
outfid = fopen(outfilename, 'a');
if outfid == -1
    error('Could not open output file');
end
if newfile
    % Output header rows
    fprintf(outfid, 'BPvsEvt\tevtIDs=\t%s\txlimalls=\t%s\tmoving_win=\t%s\n', ...
        dg_thing2str(evtIDs), dg_thing2str(xlimalls), ...
        dg_thing2str(moving_win) );
    fprintf(outfid, 'sessiondir\tlfp_FreqLim\tmaxfile');
    for k=1:length(evtIDs)
        fprintf(outfid, '\tavg');
    end
    for k=1:length(evtIDs)
        fprintf(outfid, '\tmax');
    end
    fprintf(outfid, '\n');
end
lfp_log(sprintf(...
    '\n\tStarting BPvsEvt for %s\n\tlfp_XLimAll=%s, lfp_FreqLim=%s', ...
    sessiondir, mat2str(lfp_XLimAll), mat2str(lfp_FreqLim) ));
LFPs = sort(lfp_ActiveFilenums);
for filenum = LFPs
    lfp_createWave(@lfp_bandpower, filenum, ...
        moving_win, lfp_FreqLim, ...
        'name', sprintf('%s %d-%dHz', ...
        lfp_FileNames{filenum}, lfp_FreqLim(1), lfp_FreqLim(2) ));
    for evtidx = 1:length(evtIDs)
        lfp_AlignmentRef = evtIDs{evtidx};
        lfp_XLimAll = xlimalls{evtidx};
        lfp_selectByRule('HasEvent(lfp_AlignmentRef)');
        lfp_selectByDuration('and');
        lfp_log(sprintf('Selected trials %s', ...
            dg_canonicalSeries(find(lfp_SelectedTrials))));
        [hF, data] = lfp_disp([], ...
            length(lfp_FileNames), ...
            [], 'avg');
        close(hF);
        Power(filenum, evtidx) = mean(data{1}(:,2));
    end
end
[maxRange, maxRangeIdx] = max(max(Power,[],2) - min(Power,[],2));
maxPower = Power(maxRangeIdx,:);
avgPower = mean(Power(LFPs, :), 1);
fprintf(outfid, '%s\t%s\t%g', ...
      sessiondir, dg_thing2str(lfp_FreqLim), ...
          maxRangeIdx );
fprintf(outfid, '\t%g', avgPower);
fprintf(outfid, '\t%g', maxPower);
fprintf(outfid, '\n');
fclose(outfid);
% Plot all trodes on one fig:
figure;
plot(Power(LFPs, :)','-o'); 
set(gca, 'XTick', 1:length(evtIDs)); 
xlabs={};
for k=1:length(evtIDs)
    xlabs{k} = dg_thing2str(evtIDs{k});
end
set(gca, 'XTickLabel', xlabs);
title(sprintf('BPvsEvt %s\n%d-%dHz xlimalls=%s', ...
    sessiondir, lfp_FreqLim(1), lfp_FreqLim(2), ...
    dg_thing2str(xlimalls) ));



function makeParthaPack(sessiondir)
lfp_declareGlobals;
lfp_read2('preset', sessiondir, {'Events.Nev'});
lfp_log(sprintf(...
    '\n\tStarting makeParthaPack for %s\n', ...
    sessiondir ));
numLFPs = length(lfp_ActiveFilenums);

% Load all the *.Tn and *.Tnn files for this session:
clusterfiles = {};
[animaldir, sessionname] = fileparts(sessiondir);
files = dir(animaldir);
for f = files'
    [pathstr,name,ext,versn] = fileparts(upper(f.name));
    if (length(ext) >= 2) && (ext(2) == 'T') ...
            && isequal(upper(sessionname), name(1:length(sessionname)))
        trode = str2num(ext(3:end));
        if ~isempty(trode)
            clusterfiles{end+1} = f.name;
        end
    end
end
lfp_add('preset', animaldir, ...
    clusterfiles, 'Rodent Clusters (*.Tnn, *.TTn)', false);

% Convert spikes to waves:
if length(lfp_Spikes) > 0
    for spikechannel = 1:length(lfp_Spikes)
        lfp_createWave(@lfp_spike2wave, spikechannel, ...
            'name', lfp_SpikeNames{spikechannel});
    end
end

matfilename = [lfp_SessionNames{1} 'wSpikesv6.mat'];
lfp_export(matfilename, [], '-v6');
lfp_log(['Saved ' matfilename]);




function BPvsVelo(sessiondir)
lfp_declareGlobals;

% Params for tweaking ====================================================
thinning = 100; % use 1 out of <thinning> samples for scatter plot
eventsFilename = 'Events.Nev';
bandlimitses = {[4 5], [6 10], [30 50]};
moving_win = [1 .1];
closewindows = true;
% END params for tweaking ================================================

lfp_read2('preset', sessiondir, {eventsFilename});
if isempty(lfp_XLimAll) || isempty(lfp_AlignmentRef)
    error('lfp_XLimAll and lfp_AlignmentRef must have non-empty values');
end
lfp_selectByRule('HasEvent(lfp_AlignmentRef)');
lfp_selectByDuration('and');
LFPs = sort(lfp_ActiveFilenums);
[animaldir, sessionname] = fileparts(sessiondir);
[basedir, animalID] = fileparts(animaldir);
trackerfilename = ['V' sessionname '.DAT'];
lfp_add('preset', animaldir, ...
    {trackerfilename}, 'Rodent Tracker (*.DAT)', false);
lfp_createWave(@lfp_velo, [LFPs(end) LFPs(end)+1], 'name', 'v250', 250);
vfnum = length(lfp_FileNames);
lfp_createWave(@lfp_deriv, vfnum, 'name', 'accel');
lfp_createWave(@lfp_wavesmooth, length(lfp_FileNames), 250, 'name', 'a250');
afnum = length(lfp_FileNames);
lfp_disp([ ], [afnum vfnum], [], 'fileflip', thinning, ...
    'filenames', lfp_FileNames([afnum vfnum]));
M = load('-ascii', [lfp_FileNames{vfnum} '.xls']);
velocities = reshape(M(2:end,:),[],1);
M = load('-ascii', [lfp_FileNames{afnum} '.xls']);
accelerations = reshape(M(2:end,:),[],1);
bpchannel = length(lfp_FileNames) + 1;
spreadsheetname = 'BPvsVelo.xls';
for bandidx = 1:length(bandlimitses)
    bandlimits = bandlimitses{bandidx};
    for lfpidx = 1:length(LFPs)
        lfp = LFPs(lfpidx);
        if length(lfp_FileNames) == bpchannel
            lfp_createWave(@lfp_bandpower, lfp, moving_win, bandlimits, ...
                'replace', bpchannel, 'name', 'bp');
        else
            lfp_createWave(@lfp_bandpower, lfp, moving_win, bandlimits, ...
                'name', 'bp');
        end
        lfp_disp([ ], bpchannel, [], 'fileflip', thinning, ...
            'filenames', lfp_FileNames(bpchannel));
        M = load('-ascii', [lfp_FileNames{bpchannel} '.xls']);
        bandpowers(:, lfpidx) = reshape(M(2:end,:),[],1);
    end
    avgbps = mean(bandpowers, 2);
    animalAndSession = sprintf('%s%s', animalID, sessionname);
    bandID = sprintf('%g-%gHz', bandlimits(1), bandlimits(2));
    figID = sprintf('%s %s', animalAndSession, bandID);
    titlestring = [figID ' Power vs. '];
    [hF, rV, pV] = dg_scatterplot(velocities, avgbps, [titlestring 'Velocity']);
    saveas(hF, [figID ' V.fig']);
    [hF, rA, pA] = dg_scatterplot(accelerations, avgbps, [titlestring 'Acceleration']);
    saveas(hF, [figID ' A.fig']);
    newoutfile = (exist(spreadsheetname) ~= 2);
    outfid = fopen(spreadsheetname, 'a');
    if outfid == -1
        error('Could not open output file %s for writing', ...
            spreadsheetname);
    end
    if newoutfile
        fprintf(outfid, 'session\tband\trV\tpV\trA\tpA\n');
    end
    fprintf(outfid, '%s\t%s\t%g\t%g\t%g\t%g\n', ...
        animalAndSession, bandID, rV, pV, rA, pA );
    fclose(outfid);
end
if closewindows
    close all;
end




function mycohero3(sessiondir)
% Sets lfp_SelectedTrials.
% Using the existing values of all other globals:
% Select trials that have lfp_AlignmentRef. Select trials by duration to
% maximize total analyzed time. Run coherograms for all available LFP pairs
% and save the figures, but do not print anything.

lfp_declareGlobals;

if isempty(lfp_XLimAll) || isempty(lfp_AlignmentRef)
    error('lfp_XLimAll and lfp_AlignmentRef must have non-empty values');
end

lfp_read2('preset', sessiondir, {'Events.Nev'});
lfp_log(sprintf(...
    '\n\tStarting mycohero3 for %s\n\tlfp_XLimAll=%s, lfp_FreqLim=%s, lfp_AlignmentRef=%s, lfp_CLimAll=%s', ...
    sessiondir, mat2str(lfp_XLimAll), mat2str(lfp_FreqLim), ...
    mat2str(lfp_AlignmentRef), mat2str(lfp_CLimAll) ));
numLFPs = length(lfp_ActiveFilenums);

lfp_SelectedTrials([1 2 3 4 5 6 10 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40])=false;

%lfp_SelectedTrials([5 7 8 10 19 22 23 24 25 27 32 39 40])=false;
%lfp_selectByRule('HasEvent([17 18])');
%lfp_BadTrials = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 20 29 37 41 42 43 44 45 47 52];
%lfp_BadTrials = [1 2 3 4 5 6 7 8 10 11 13 14 17 37 41 45 47];
%lfp_BadTrials=union(union([6 7 8 12 13 22 25 29 32 36 37],[12 13 28]),[20 28 37]);
% lfp_selectByDuration('and');
lfp_log(sprintf('Selected trials %s', dg_canonicalSeries(find(lfp_SelectedTrials))));
if ~any(lfp_SelectedTrials)
    return
end

% Params for tweaking ====================================================
moving_win = [1 .1];
nw = 1.8;
k = 3;
padding = 2;
% END params for tweaking ================================================

[animaldir, sessionname] = fileparts(sessiondir);
[basedir, animalID] = fileparts(animaldir);
filenumpairs = [];
LFPs = sort(lfp_ActiveFilenums);
for lfpA = LFPs(1:end-1)
    for lfpB = LFPs(find(LFPs > lfpA))
        lfp_log(sprintf('Coherogram: %s vs %s', ...
            lfp_FileNames{lfpA}, ...
            lfp_FileNames{lfpB} ));
        hF = lfp_spec('coh', [], [lfpA lfpB], moving_win, ...
            'rmdc', 'nw', nw, 'k', k, 'logsig', ...
            'avg', 'pad', padding);
        saveas(hF, sprintf('%s%s_%s%s_%d.fig', ...
            animalID, sessionname, ...
            lfp_FileNames{lfpA}, lfp_FileNames{lfpB}, ...
            lfp_AlignmentRef(1) ));
        close(hF);
    end
end




function makeParthaPack2(sessiondir)
lfp_declareGlobals;
lfp_read2('preset', sessiondir, {'Events.Nev'});
lfp_log(sprintf(...
    '\n\tStarting makeParthaPack for %s\n', ...
    sessiondir ));
numLFPs = length(lfp_ActiveFilenums);
if ismember(lfp_AlignmentRef, [23 14])
    lfp_XLimAll = [-1 5];
elseif isequal(lfp_AlignmentRef, [15 16]) ...
        || isequal(lfp_AlignmentRef, [17 18])
    lfp_XLimAll = [-4 2];
end
lfp_selectByRule('HasEvent(lfp_AlignmentRef)');
lfp_selectByDuration('and');
lfp_log(sprintf(...
    '\n\tSelected trials %s\n', ...
    dg_canonicalSeries(find(lfp_SelectedTrials)) ));
index = lfp_findLastPointByValue([], [], 2047, -2048, 3);
lfp_XLimAll(1) = ceil(4*index*lfp_SamplePeriod)/4;
lfp_excludeByValue([], [], 2047, -2048);

% Load all the *.Tn and *.Tnn files for this session:
clusterfiles = {};
[animaldir, sessionname] = fileparts(sessiondir);
files = dir(animaldir);
for f = files'
    [pathstr,name,ext,versn] = fileparts(upper(f.name));
    if (length(ext) >= 2) && (ext(2) == 'T') ...
            && isequal(upper(sessionname), name(1:length(sessionname)))
        trode = str2num(ext(3:end));
        if ~isempty(trode)
            clusterfiles{end+1} = f.name;
        end
    end
end
lfp_add('preset', animaldir, ...
    clusterfiles, 'Rodent Clusters (*.Tnn, *.TTn)', false);

% Convert spikes to waves:
if length(lfp_Spikes) > 0
    for spikechannel = 1:length(lfp_Spikes)
        lfp_createWave(@lfp_spike2wave, spikechannel, ...
            'name', lfp_SpikeNames{spikechannel});
    end
end

% De-select any trials that end up empty after applying lfp_XLimAll:
for trial = lfp_enabledTrials(1:size(lfp_TrialIndex,1))
    [samplerange, refpoint] = lfp_getSampleRange(trial);
    if isempty(samplerange)
        lfp_SelectedTrials(trial) = false;
    end
end

matfilename = sprintf('%swSp_e%d_v6.mat', ...
    lfp_SessionNames{1}, lfp_AlignmentRef(1) );
if sum(lfp_SelectedTrials) == 0
    emptymsg = sprintf(...
        'Skipping session %s because it contains no usable data', ...
        lfp_SessionNames{1} );
    warning(emptymsg);
    lfp_log(emptymsg);
else
    lfp_export(matfilename, [], '-v6', 'selectedTrials', 'useXLim');
    lfp_log(['Saved ' matfilename]);
end



