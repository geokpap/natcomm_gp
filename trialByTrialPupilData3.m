function trialByTrialPupilData3(sessiondir)
% Per Georgios' email of 1/17/24, 7:14 PM.  Pupil diameter data only.
%INPUT
% sessiondir
%OUTPUT
% All output goes to a "trialByTrialPupilData.mat" file in the directory in
% the <outroot> tree that corresponds to <sessiondir> as in
% 'cleanEyeLinkSession'.  The output file contains fields 'choicedata' and
% 'forceddata'.
%NOTES
% Intermediate storage is in <outroot>, and aggregation and addition of
% animal and session IDs must be done in a second pass.  First we check to
% see if the cleaned, downsampled data exist in <outroot>, and if not we
% find and copy it from wherever it is in <altroots>.  If it is not there
% either, then we run 'cleanEyeLinkSession2'.
%   Sessiondirs that have 'fixed' subdir are renamed on output by appending
% 'fixed' to the sessionID and eliminating the subdir.
%   AnimalIDs are standardized on output to {Debbie|Prez}.
%   Only correctly completed trials (i.e. with outcome events) are
% included.
%   Wed Apr  3 15:02:38 CST 2024: Only cleans csc126 if the 
% 'csc126clean_down32.mat' file does not exist.
%   Thu Apr  4 14:23:00 CST 2024: Added call to 'getGPSessionType'.  Added
% code to skip any ApAv (or ApAv_em) sessions that already have a
% "trialByTrialPupilData.mat" file.  Modified code to raise a warning and
% replace the output file for other session types.
%   Sep 6-8 2024: assorted tweaks.

global lfp_Samples lfp_Events lfp_TrialIndex lfp_TrialParams %#ok<GVMIS>

outroot = '/annex2/analysis/dgibson/pupilApAv_rev_only'; % as in 'cleanEyeLinkSession'
altroots = {
    '/annex2/analysis/dgibson/eyelink'
    '/annex2/analysis/dgibson/eyelink2'
    };
outfilename = 'trialByTrialPupilData.mat';
blinklev = -0.05;
prezIDs = {'prez' 'Prez'};
debIDs = {'debbie' 'LittleDebbie'};

sestype = getGPSessionType(sessiondir);
lfp_changeSetup(sestype);
sesIDsuffix = '';
pathparts = strsplit(sessiondir, '/');
if isequal(pathparts{end}, 'fixed')
    pathparts(end) = [];
    sesIDsuffix = 'fixed';
end
if isempty(pathparts{1})
    pathparts(1) = [];
end
% <sessiondir> formats to handle:
%   /annex?/georgios/animalID/ApAv/sessionID
%   /annex?/georgios/animalID/{Recordings|Behavior}/ApAv/sessionID
%   /annex?/georgios/animalID/...
if ismember(pathparts{3}, [prezIDs debIDs])
    if ismember(pathparts{3}, prezIDs)
        animalID = 'Prez';
    else
        animalID = 'Debbie';
    end
else
    error('trialByTrialPupilData3:pathname', ...
        'Unknown pathname format: %s', sessiondir);
end

outdir = fullfile(outroot, animalID, pathparts{4:end});
fprintf('outdir = %s\n', outdir);
if exist(outdir, 'dir')
    if exist(fullfile(outdir, outfilename), 'file')
        if isequal(sestype, 'georgios_ApAv')
            % do nothing
            return
        else
            warning('trialByTrialPupilData:outfile', ...
                'Replacing output file.');
        end
    end
else
    dg_mkdir_p(outdir);
end

% if necessary, copy 'csc126clean_down32.mat' with evts file
newoutdir = [outdir sesIDsuffix];
if ~exist( fullfile(newoutdir, 'csc126clean_down32.mat'), ...
        'file' )
    srcdir = ''; % where 'csc126clean_down32.mat' actually is
    success = false;
    for altidx = 1:length(altroots)
        srcdir = fullfile(altroots{altidx}, animalID, pathparts{4:end});
        if exist(fullfile([srcdir sesIDsuffix], ...
                'csc126clean_down32.mat'), 'file')
            success = true;
            break
        end
        if exist(fullfile(srcdir, ...
                'csc126clean_down32.mat'), 'file')
            success = true;
            break
        end
    end
    if ~success
        % There is no cleaned downsampled file!
        cleanEyeLinkSession2(sessiondir);
        if ~isempty(sesIDsuffix)
            % Append <sesIDsuffix> to <outdir> and move the files from the
            % original <sesIDsuffix> subdirectory to the new <outdir>.
            [status,msg] = movefile(outdir, newoutdir);
            if ~status
                error('trialByTrialPupilData3:outdir', ...
                    'Could not move output from cleanEyeLinkSession in %s\n%s', ...
                    outdir, msg);
            end
            [status,msg] = movefile(fullfile(newoutdir, sesIDsuffix, '*'), ...
                newoutdir);
            if ~status
                error('trialByTrialPupilData3:sesIDsuffix', ...
                    'Could not move output from sesIDsuffix dir %s\n%s', ...
                    fullfile(newoutdir, sesIDsuffix, '*'), msg);
            end
            [status,msg] = rmdir(fullfile(newoutdir, sesIDsuffix));
            if ~status
                error('trialByTrialPupilData3:sesIDsuffixdir', ...
                    'Could not remove sesIDsuffix dir %s\n%s', ...
                    fullfile(newoutdir, sesIDsuffix), msg);
            end
        end
    end
    if ~exist(newoutdir, 'dir')
        mkdir(newoutdir);
    end
    dg_copyfile( fullfile( [srcdir sesIDsuffix], ...
        'csc126clean_down32.mat' ), ...
        fullfile(newoutdir, 'csc126clean_down32.mat') );
    if exist(fullfile([srcdir sesIDsuffix], 'events.mat'), 'file')
        dg_copyfile( fullfile([srcdir sesIDsuffix], 'events.mat'), ...
            fullfile(newoutdir, 'events.mat') );
    elseif exist(fullfile([srcdir sesIDsuffix], 'Events.mat'), 'file')
        dg_copyfile( fullfile([srcdir sesIDsuffix], 'Events.mat'), ...
            fullfile(newoutdir, 'events.mat') );
    else
        error('trialByTrialPupilData3:oops', 'Can''t copy!');
    end
end
outdir = newoutdir;

% We now have cleaned pupil diameter output to use.  Replace blink-level
% values with NaN after reading.
lfp_getEvtIDs;
filenames = {'events.mat' 'csc126clean_down32.mat'};
dirfiles = dir(outdir);
dirfilenames = {dirfiles(~[dirfiles.isdir]).name};
infiles = cell(0, 2);
for fidx = 1:length(filenames)
    if ismember(filenames{fidx}, lower(dirfilenames))
        infiles{fidx} = dirfilenames{ ...
            ismember(lower(dirfilenames), filenames{fidx}) };
    else
        error('trialByTrialPupilData3:missingfile', ...
            'There is no %s in %s', filenames{fidx}, outdir);
    end
end
lfp_read2('preset', outdir, infiles);
lfp_Samples{1}(lfp_Samples{1}(:) < blinklev) = NaN;

% Compute average pupil diameter for the time between the cue onset and
% right at the first outcome delivery i.e. the onset of the airpuff for
% choosing Approach, or the onset of the small reward for choosing
% avoidance.
%   1.  Choice trials:
lfp_selectByRule( ...
    'HasEvent([CROSS_LEFT CROSS_RIGHT]) && HasEvent([airpuffOnAp rewardOnAv])');
numchoice = length(lfp_enabledTrials);
choicedata = NaN(numchoice, 6);
rownum = 1;
for tnum = lfp_enabledTrials
    evts = lfp_Events(lfp_TrialIndex(tnum, 1):lfp_TrialIndex(tnum, 2), :);
    cuesamp = lfp_time2index(evts(find(evts(:,2) == choiceCueOn, 1)  ));
    outcomesamp = lfp_time2index( ...
        evts(find(ismember(evts(:,2), [airpuffOnAp rewardOnAv]), 1)) );
    magred = lfp_TrialParams{tnum}(1);
    magyel = lfp_TrialParams{tnum}(2);
    isAp = any(evts(:,2) == airpuffOnAp);
    avgdiam = nanmean(lfp_Samples{1}(cuesamp:outcomesamp)); %#ok<NANMEAN>
    meddiam = nanmedian(lfp_Samples{1}(cuesamp:outcomesamp)); %#ok<NANMEDIAN>
    trialTS = lfp_Events(lfp_TrialIndex(tnum, 1), 1);
    choicedata(rownum, :) = [magred, magyel, isAp, avgdiam, ...
        meddiam, trialTS];
    rownum = rownum + 1;
end
%   2.  Forced (Pavlovian) trials:
lfp_selectByRule( ...
    'HasEvent(FORCED) && HasEvent([AVE_ON_FORCED RWD_ON_FORCED])');
numforced = length(lfp_enabledTrials);
forceddata = NaN(numforced, 5);
rownum = 1;
for tnum = lfp_enabledTrials
    evts = lfp_Events(lfp_TrialIndex(tnum, 1):lfp_TrialIndex(tnum, 2), :);
    cuesamp = lfp_time2index(evts( ...
        find(ismember(evts(:,2), [forcedRedOn forcedYelOn]), 1) ));
    outcomesamp = lfp_time2index(evts( ...
        find(ismember(evts(:,2), [AVE_ON_FORCED RWD_ON_FORCED]), 1) ));
    if lfp_TrialParams{tnum}(1) ~= 0 && lfp_TrialParams{tnum}(2) == 0
        isrwd = true;
        mag = lfp_TrialParams{tnum}(1);
    elseif lfp_TrialParams{tnum}(2) ~= 0 && lfp_TrialParams{tnum}(1) == 0
        isrwd = false;
        mag = lfp_TrialParams{tnum}(2);
    else
        error('trialByTrialPupilData3:parms', ...
            'Invalid param pair for forced trial, trialnum %d.', tnum);
    end
    avgdiam = nanmean(lfp_Samples{1}(cuesamp:outcomesamp)); %#ok<NANMEAN>
    meddiam = nanmedian(lfp_Samples{1}(cuesamp:outcomesamp)); %#ok<NANMEDIAN>
    trialTS = lfp_Events(lfp_TrialIndex(tnum, 1), 1);
    forceddata(rownum, :) = [mag, isrwd, avgdiam, meddiam, trialTS];
    rownum = rownum + 1;
end
save(fullfile(outdir, outfilename), 'choicedata', 'forceddata');
fprintf('Done saving %s\n', fullfile(outdir, outfilename));
