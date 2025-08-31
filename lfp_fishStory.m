function lfp_fishStory(tackle, sessiondir, varargin)
% A quick analysis of putative striosomal units following the Friedman
% protocol.  Inserts lfp_NominalTrialStart events 10 us before each
% stimulation or baseline event, and lfp_NominalTrialEnd events 1.99 s
% after each stimulation or baseline event.  If stimulations are 2 s apart,
% this results in an ITI of 9.99 ms.  The baseline period is defined to be
% the period from the first recorded event until 10 us before the first
% stimulation event.  Files in <sessiondir> are assumed to be CSC files if
% and only if the filename matches "csc_*".
%INPUTS
% tackle:  a structure containing everything you need for fishing, in the
%   following fields: 
%   artifactTimeout: the number of microseconds over which to linearly
%       interpolate LFPs after each stimulation event.
%   baselineID: the TTL event code for markers in baseline period, without
%       strobe. 
%   baselineT: the time period between markers during the baseline period,
%       in seconds.
%   baselineDelay: the time period between the first recorded event and the
%       first marker in the baseline period, in seconds.
%   blkIDs: a vector of unstrobed TTL event codes for stimulation in each
%       block, where a series of the same evt code defines the block.
%       These must be specified in block numeric order, i.e. the first
%       value defines block 1, the second block 2, etc.
%   downsampledRate: the downsampling factor chosen automatically for
%        dg_downsampleAndConvert is the largest integer that yields a
%        downsampled sample rate of at least downsampledRate Hz.  It is
%        assumed that all CSC files are sampled at the same rate.
%   evfilename: name of events file relative to <sessiondir>.
%   localavgrefs: a cell vector of strings or cell arrays of strings where
%       each element is a legitimate value for the <CSCfilespec> argument
%       to dg_localAvgRef that specifies the files to average together and
%       save as a "local average reference" channel.  The results are saved
%       to files named localavgref<N>, where <N> =
%       1:length(tackle.localavgrefs).
%   spikefiles: a value suitable for submission to dg_spikes2nev, or empty
%       if the events file is already perfect just the way it is.
% sessiondir: pathname to directory containing data to analyze
%OUTPUTS
% ...are all files saved to <sessiondir>.  If lfp_fishStory.nev already
%   exists, it gets overwritten.  If other files already exist, creation
%   of those files gets skipped.
% 'lfp_fishStory.nev': the merged events file, containing
%   lfp_NominalTrialStart events at 10 us before each stimulation/baseline,
%   and lfp_NominalTrialEnd events at 20 us before the second through last
%   stimulation/baselines, and a final lfp_NominalTrialEnd placed to make
%   the last trial have a duration equal to the median of its predecessors.
% 'baseline': fragment directory
% 'block<n>': fragment directory from stimulation site # <n>
%OPTIONAL FIELDS IN <tackle>
% 'startevt' - search for an event with TTL value 0 and Event String value
%   that case-insensitively matches <tackle.startevt> to mark start of
%   baseline period instead of the first recorded event.
% 'baselinedur' - analyze only the last <tackle.baselinedur> seconds of the
%   baseline period.  If there are not <tackle.baselinedur> seconds of
%   recording available in the baseline period, a warning is issued and the
%   entire baseline period is used.


%$Rev: 262 $
%$Date: 2012-02-09 14:00:54 -0500 (Thu, 09 Feb 2012) $
%$Author: dgibson $

lfp_declareGlobals;

% Merge and create event markers as needed:
evpath = fullfile(sessiondir, tackle.evfilename);
mergedeventfname = 'lfp_fishStory.nev';
mergedeventpath = fullfile(sessiondir, mergedeventfname);
if isempty(tackle.spikefiles)
    dg_copyfile(evpath, mergedeventpath);
else
    dg_spikes2nev(mergedeventpath, tackle.spikefiles, {evpath});
end
[evTS, TTL, ES, Hdr] = dg_readEvents(mergedeventpath); 

if isempty(evTS)
    error('lfp_fishStory:noevt', ...
        'There are no events in the session.');
end
TTL(TTL<0) = TTL(TTL<0) + 2^15;
if ismember(TTL(1), tackle.blkIDs)
    firstevt = [];
else
    firstevt = [evTS(1) TTL(1)];
end
evtidx = cell(size(tackle.blkIDs));
for blkidx = 1:length(tackle.blkIDs)
    evtidx{blkidx} = find(TTL==tackle.blkIDs(blkidx));
end
BLstartevtidx = 1;
if isfield(tackle, 'startevt')
    evt0idx = find(TTL==0);
    if isempty(evt0idx)
        warning('lfp_fishStory:evt0idx', ...
            'There is no TTL=0, using first event for startevt.');
    else
        matchidx = find(cellfun( @strcmpi, ES, ...
            repmat({tackle.startevt}, size(ES)) ));
        if isempty(matchidx)
            warning('lfp_fishStory:matchidx', ...
                'There is no "%s" event, using first event for startevt.', ...
                tackle.startevt );
        else
            BLstartevtidx = evt0idx(matchidx(1));
        end
    end
end
% Event file timestamps are in microseconds, so convert tackle timing
% params to microseconds:
BLendTS = evTS(evtidx{1}(1)) - tackle.baselineT * 1e6;
if ~isfield(tackle, 'baselinedur') || ...
        tackle.baselinedur * 1e6 > evTS(evtidx{1}(1)) - evTS(BLstartevtidx)
    if isfield(tackle, 'baselinedur')
        warning( 'lfp_fishStory:baselinedur', ...
            'There are only %.6f seconds of baseline time available, using all.', ...
            (evTS(evtidx{1}(1)) - evTS(BLstartevtidx))*1e-6 );
    end
    BLstartTS = evTS(BLstartevtidx) + tackle.baselineDelay * 1e6;
else
    BLstartTS = evTS(evtidx{1}(1)) - tackle.baselinedur * 1e6;
end
numBL = floor((BLendTS - BLstartTS) / (tackle.baselineT * 1e6));
BLTS = BLstartTS + (0 : numBL - 1) * (tackle.baselineT * 1e6);
newTS = BLTS;
newTTL = repmat(tackle.baselineID, size(BLTS));
for blkidx = 1:length(tackle.blkIDs)
    if isempty(evtidx{blkidx})
        warning('lfp_fishStory:nostim', ...
            'There are no stimulation events with ID %d', ...
            tackle.blkIDs(blkidx));
    else
        newTS = [newTS evTS(evtidx{blkidx})]; %#ok<AGROW>
        newTTL = [newTTL TTL(evtidx{blkidx})]; %#ok<AGROW>
    end
end
% At this point, newTS and newTTL contain only stimulation and baseline events. 
lfp_loadSetup;
newTS = [ newTS newTS-10 newTS(2:end)-20 newTS(end)+median(diff(newTS)) ];
newTTL = [ newTTL lfp_NominalTrialStart(1)*ones(size(newTTL)) ...
    lfp_NominalTrialEnd(1)*ones(size(newTTL)) ];
[newTS, idx] = sort(newTS);
newTTL = newTTL(idx);
% This must be done last so the "first event" doesn't become a trial unto
% itself:
if ~isempty(firstevt)
    newTS = [firstevt(1) newTS];
    newTTL = [firstevt(2) newTTL];
end

% Write the final event file:
Hdr{2,1} = sprintf( ...
    '## File Name: %s ', mergedeventpath);
Hdr{4,1} = sprintf( ...
    '## Time Closed: (mm/dd/yyyy): %s  At Time (HH:MM:SS): %s ', ...
    datestr(now, 23), datestr(now, 13));
dg_writeNlxEvents(mergedeventpath, newTS, newTTL - 2^15, {}, Hdr);
fprintf('Wrote %s\n', mergedeventpath);

% We now have the final set of events in newTS, newTTL.  This is the time
% to remove the stimulation artifacts and downsample.

CSCfiles = dir(fullfile(sessiondir, 'csc_*'));
dest = fullfile(sessiondir, 'zapped');
if exist(dest) %#ok<*EXIST>
    if exist(dest) ~= 7
        error('There is a file named "%s"', dest);
    end
else
    mkdir(dest);
end
if isempty(CSCfiles)
    error('lfp_fishStory:csc', ...
        'There are no CSC files in %s', sessiondir);
end
sampleperiod = [];
for fileidx = 1:length(CSCfiles)
    % Throughout this loop, we assume that tackle.blkIDs lists the event
    % IDs in chronological order of occurence, i.e. that block 1 is truly
    % the first, block 2 is the second, etc.  We also assume that 5 s of
    % elapsed time is a reasonable range over which to search for the next
    % stimulation event.
    zapfilepath = fullfile(dest, CSCfiles(fileidx).name);
    if exist(zapfilepath, 'file')
        warning('lfp_fishStory:zapfile', ...
            'File "%s" already exists, skipping', zapfilepath);
        continue
    end
    [cscTS, Samples, Hdr] = dg_readCSC( ...
        fullfile(sessiondir, CSCfiles(fileidx).name) );
    sampleperiod = median(diff(cscTS)) / size(Samples,1);
    numIntrp = round(tackle.artifactTimeout / sampleperiod);
    allTS = repmat(cscTS, size(Samples,1), 1) + sampleperiod * ...
        repmat((0:size(Samples,1) - 1)', 1, size(Samples, 2));
    laststimsample = 1;
    for blkidx = 1:length(tackle.blkIDs)
        stimevtidx = find(newTTL==tackle.blkIDs(blkidx));
        for stimnum = 1:length(stimevtidx)
            startsample = laststimsample;
            endsample = min( numel(Samples), ...
                startsample + floor(5e6/sampleperiod) );
            index = dg_binsearch(allTS(:), newTS(stimevtidx(stimnum)), ...
                startsample, endsample);
            if index < numel(allTS)+1
                if allTS(index) - newTS(stimevtidx(stimnum)) > sampleperiod/2
                    index = index - 1;
                end
                Samples(index + (0:numIntrp)) = linspace(Samples(index), ...
                    Samples(index+numIntrp), numIntrp+1);
            else
                break
            end
            laststimsample = index;
        end
    end
    dg_writeCSC(zapfilepath,...
        cscTS, Samples, Hdr);
    fprintf('Wrote %s\n', fullfile(dest, CSCfiles(fileidx).name));
end
dg_copyfile(mergedeventpath, dest);

% Now downsample:
src = fullfile(sessiondir, 'zapped');
dsdest = fullfile(sessiondir, 'downsampled');
if exist(dsdest)
    if exist(dsdest) ~= 7
        error('There is a file named "%s"', dsdest);
    end
else
    mkdir(dsdest);
end
if exist(dsdest) && exist(dsdest) ~= 7
    error('There is a file named "%s"', dsdest);
end
fprintf('Downsampling...');
if isempty(sampleperiod)
    [cscTS, Samples] = dg_readCSC( ...
        fullfile(sessiondir, CSCfiles(1).name) );
    sampleperiod = median(diff(cscTS)) / size(Samples,1);
end
ratefactor = floor(1 / (sampleperiod * 1e-6 * tackle.downsampledRate) );
dg_downsampleAndConvert(src, dsdest, ratefactor, ...
    'reconcile', 'all', 'skip');
fprintf(' done\n');

% Compute local average refs:
for refidx = 1:length(tackle.localavgrefs)
    CSCfilespec = tackle.localavgrefs{refidx};
    outfilepath = fullfile(dsdest, sprintf( ...
        'csc_localavgref%d.mat', refidx ));
    dg_localAvgRef(dsdest, CSCfilespec, outfilepath);
end

save(fullfile(dsdest, 'lfp_fishStory_tackle.mat'), 'tackle');

% Fragment the session:
% (Do we still need this with downsampled files?)
disp('Starting fragmentation');
rules = {'baseline'	sprintf('HasEvent(%d)', tackle.baselineID)};
for blkidx = 1:length(tackle.blkIDs)
    rules(end+1,:) = {sprintf('block%d', blkidx) ...
        sprintf('HasEvent(%d)', tackle.blkIDs(blkidx))}; %#ok<AGROW>
end
[p, name] = fileparts(mergedeventfname); %#ok<ASGLU>
lfp_fragmentFiles(rules, 'preset', dsdest, {[name '.mat']});
fprintf('Fragmented %s\n', sessiondir);

