function bulkProcess_bill(processnum, dirlist, varargin)
%BULKPROCESS_BILL performs offline lfp_lib processing on a list of
%directories.
% bulkProcess_joey(processnum, dirlist)
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
%       only processes units listed in lookup table

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
    error('lfp_XLimAll and lfp_AlignmentRef must have values for ''mycoherograms''');
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
    error('lfp_XLimAll and lfp_AlignmentRef must have values');
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
    error('lfp_XLimAll and lfp_AlignmentRef must have values');
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



function mycohero2(sessiondir)
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
    error('lfp_XLimAll and lfp_AlignmentRef must have values');
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
        lfp_spec('coh', [], [lfp_ActiveFilenums(lfpnum) lfp_ActiveFilenums(end)], moving_win, ...
            'rmdc', 'nw', nw, 'k', k, 'logsig', ...
            'avg', 'minspikes', minspikes);
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
    error('lfp_XLimAll and lfp_AlignmentRef must have values');
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
    'C:\B08\Acq24' {'ACQ24-T2C4','ACQ24-T6C3'}
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
    % remove EP from spike wave:
    [hF, data] = lfp_disp([], ...
        [wavechannel], ...
        [], 'avg');
    close(hF);
    lfp_createWave(@lfp_rmEP, wavechannel, ...
        data{1}(1,1), data{1}(:,2), ...
        'name', [lfp_FileNames{wavechannel} 'rmEP']);
    newwavechannel = length(lfp_FileNames);
    for lfpnum = 1:numLFPs
        % remove EP from LFP:
        [hF, data] = lfp_disp([], ...
            [lfp_ActiveFilenums(lfpnum)], ...
            [], 'avg');
        close(hF);
        lfp_createWave(@lfp_rmEP, lfp_ActiveFilenums(lfpnum), ...
            data{1}(1,1), data{1}(:,2), ...
            'name', [lfp_FileNames{lfp_ActiveFilenums(lfpnum)} 'rmEP']);
        newlfpchannel = length(lfp_FileNames);
        % log a message and compute coherence
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

