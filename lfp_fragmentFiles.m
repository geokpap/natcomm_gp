function lfp_fragmentFiles(rules, varargin)
%lfp_fragmentFiles(rules)
%lfp_fragmentFiles(rules, 'preset', presetDir, presetFiles)
% <rules> is a two-column cell array.  The first column contains a name for
% each rule.  Each element in the second column is either a matrix of
% parameter numbers and values for submission to lfp_trialHasParams, or a
% string containing an arbitrary Matlab expression that evaluates to true
% for the trials to be included (see lfp_selectByRule for details).
%OPTIONS
%   'destdir', destdir - creates output files in a different directory tree
%       from lfp_DataDir.  The value of <destdir> should normally be an
%       animal directory.  The last element in the path lfp_DataDir is
%       taken to be a session directory name, and the rule-named
%       subdirectories are created in a directory of the same name which is
%       in turn a subdirectory of <destdir>.  In other words, if you say
%       "'destdir', 'c:\mugwump'", and lfp_DataDir (or presetDir when using
%       the 'preset' option) is 'd:\mydata\monkeys\mugwump\20090101', then
%       the fragment subdirectories will be created in
%       'c:\mugwump\20090101'.
%   'preset', presetDir, presetFiles -
%       This works as for lfp_read2, i.e. the GUI is bypassed, <presetDir>
%       specifies the data directory, <presetFiles> specifies the files to
%       create fragments from, the first of which is the events file.  If
%       <presetFiles> contains only the events file name, then CSC files
%       are identified automatically as when using lfp_read2 with
%       lfp_UseFileSelect=false (see comments on CSCFileRegexp, CSCFileExt
%       in lfp_preamble).
%   'noNEV' - skips fragmentation of NEV files (except for the events file
%       specified, if it ends in .NEV).

% The general strategy is to select the trials that satisfy one of the
% <rules>, and then determine from lfp_TrialIndex what time intervals are
% needed to cover those trials.   It is up to the user to ensure (or not)
% that the rules adequately cover the input data.  See header comments for
% function findTimeIntervals in this file for further details on which time
% intervals and events will or will not be included.  Time intervals are
% expanded to fill frame boundaries when fragmenting CSC files.  
%
% A log file named lfp_fragmentFiles.log is maintained.  Type 'which
% lfp_fragmentFiles.log' to find out what directory it's in.  If there is
% no such file, a new one will be started in your current Working
% Directory.
%
% Each output directory contains a file named 'lfp_FragmentFiles.mat' that
% contains the session name and the list of trials in the current fragment.
%
% A new subdirectory with the name of the rule is
% created in lfp_DataDir to contain the output files.  All *.NEV and *.NCS
% files are read as input. The output files have the same base names as the
% corresponding input files, but they are Matlab data files with the .MAT
% extension. Only the events that are within the originally
% determined time intervals (before expanding to fill frames) are written
% to the output.
%
% In CSC files,
% dg_Nlx2Mat_Samples is frame formatted in points X frames format as
% returned by Nlx2MatCSC, and dg_Nlx2Mat_Timestamps is a row vector as
% returned by Nlx2MatCSC.  In event files (not to be confused with .EVTSAV
% files), dg_Nlx2Mat_EventStrings is a char matrix with one event strings
% per row, as obtained by applying cell2mat to the result returned by
% Nlx2MatEV, and (consistently with the Nlx functions)
% dg_Nlx2Mat_Timestamps is a row vector as for CSC files.  See lfp_save for
% specs of .EVTSAV files.  Note that the events.MAT files that result from
% fragmenting an events.NEV file contain only timestamps and event strings,
% so it is essential that the event strings be in the "standard" (i.e.
% obsolete as of Cheetah 5) Neuralynx format.  This problem can be
% circumvented by saving an .EVTSAV file before fragmenting, or by saving a
% .NEV file with standard event strings produced by calling
% dg_writeNlxEvents with an empty event string array.
%
% Except in the case of .EVTSAV files, timestamps are used throughout
% exactly as they were read, so time units in the output files are the same
% as in the original files. For Neuralynx files that would be microseconds.
% In the case of .EVTSAV files, the timestamps are multiplied by 1e6 on
% reading so that they will be compatible with the timestamps in Neuralynx
% files, and then they are divided by 1e6 on output so that they will be in
% spec for .EVTSAV files. Also, samples in CSC files are formatted with the
% same number of samples per frame as in the input file.
%
% lfp_NewTrialTarget invokes the use of ALREADY EXISTING lfp_NewTrialStart
% and lfp_NewTrialEnd events, but if they will not be created if they don't
% exist (unlike lfp_read2 and lfp_add).

%$Rev: 418 $
%$Date: 2022-08-06 20:48:17 -0400 (Sat, 06 Aug 2022) $
%$Author: dgibson $

global lfp_fragmentFiles_save_events lfp_fragmentFiles_save_params
% Initialize those globals to be empty (clear does NOT do this!):
lfp_fragmentFiles_save_events = [];
lfp_fragmentFiles_save_params = [];

argnum = 1;
destdir = '';
fragNEV = true;
presetflag = false;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'destdir'
            argnum = argnum + 1;
            destdir = varargin{argnum};
        case 'preset'
            presetflag = true;
            presetDir = varargin{argnum+1};
            presetFiles = varargin{argnum+2};
            argnum = argnum + 2;
        case 'noNEV'
            fragNEV = false;
        otherwise
            error('lfp_fragmentFiles:badoption', ...
                'The option "%s" is not recognized.', ...
                dg_thing2str(varargin{argnum} ) );
    end
    argnum = argnum + 1;
end


lfp_preamble;
if newdatadir == 0
    return
end

lfp_SamplesPerFrame = 512;

if presetflag
    % EVFilename2 comes from setup, but lfp_readEvents will raise an error
    % in 'preset' mode if it can't find EVFilename1.
    EVFilename1 = presetFiles{1};
    % Construct CSCfilelist:
    if length(presetFiles) > 1
        CSCfilelist = presetFiles(2:end);
    else
        CSCDataFiles = dir(lfp_DataDir);
        CSCfilelist = {};
        for entrynum = 1:length(CSCDataFiles)
            file = CSCDataFiles(entrynum);
            if file.bytes == 16384
                % empty file; unceremoniously ignore
                continue
            end
            filenameCSC = file.name;
            pathnameCSC = fullfile(lfp_DataDir, filenameCSC);
            % This is where we decide which files to open:
            if regexpi(filenameCSC, [CSCFileRegexp '\.' CSCFileExt])
                CSCfilelist{end+1, 1} = filenameCSC;
            end
        end
    end
    newdatadir =  presetDir;
    read_mode = 'preset';
else
    % NOTE:  CSCFileRegExp is not used in this case.
    % Override filename values from the setup file.  These values only
    % matter if ~lfp_UseFileSelect.
    EVFilename1 = 'Events.nev';
    EVFilename2 = 'Events.dat';
    read_mode = 'nlx';
end

[lfp_TrialParams, lfp_Events, fnameEV, badtrials] = ...
    lfp_readEvents( read_mode, EVFilename1, EVFilename2, ...
    false, false, usestrobe );
if ~isempty(badtrials)
    error('lfp_fragmentFiles:badtrials', ...
        'There are bad trials reported by ''lfp_readEvents''.');
end
if isempty(lfp_TrialParams) && isempty(lfp_Events) && isempty(fnameEV)
    return;
end
[pathstr,name,ext] = fileparts(fnameEV);
fnameEV = [name ext];

% Find a CSC file, read the timestamps, and calculate lfp_SamplePeriod,
% leaving all time measurements in raw Neuralynx units:
pathname = fullfile(lfp_DataDir, CSCfilelist{1});
[pathstr,name,ext] = fileparts(pathname);
switch (upper(ext))
    case '.MAT'
        load('-mat', pathname);
        lfp_TimeStamps = dg_Nlx2Mat_Timestamps;
        clear dg_Nlx2Mat_Timestamps;
        clear dg_Nlx2Mat_Samples;
    otherwise
        lfp_TimeStamps = dg_readCSC(pathname);
end
lfp_SamplePeriod = median(diff(lfp_TimeStamps)) ...
    / lfp_SamplesPerFrame;

[origtrialnums, sessionname] = lfp_getSessionInfo;

% Use new start and end events if specified in setup:
if ~isempty(lfp_NewTrialTarget)
    lfp_NominalTrialStart = lfp_NewTrialStart;
    lfp_NominalTrialEnd = lfp_NewTrialEnd;
end

lfp_TrialIndex = lfp_createTrialIndex(lfp_Events);

lfp_SelectedTrials = repmat(true, 1, size(lfp_TrialIndex,1));

for rulenum = 1:size(rules,1)
    if exist('lfp_fragmentFiles.log') == 2
        logfname = which('lfp_fragmentFiles.log');
    else
        logfname = 'lfp_fragmentFiles.log';
    end
    logfid = fopen(logfname, 'a');
    if logfid < 0
        warning('lfp_fragmentFiles:nologfile', ...
            'Could not open log file.' );
    else
        fprintf(logfid, '%s rulenum %d %s %s\n', ...
            datestr(now, 0), rulenum, ...
            dg_thing2str(rules{rulenum,1}), dg_thing2str(rules{rulenum,2}) );
        fclose(logfid);
    end
    try lfp_selectByRule(rules{rulenum,2}, rules{rulenum,1});
    catch
        disp(sprintf('Rule number %d, ''%s'':', rulenum, rules{rulenum,1}));
        disp(lasterr);
        continue
    end
    selectedtime = findTimeIntervals;
    createFragments(rules{rulenum,1}, selectedtime, CSCfilelist, fnameEV, ...
        origtrialnums, sessionname, fragNEV, destdir);
end

clear lfp_fragmentFiles_save_events lfp_fragmentFiles_save_params
disp('lfp_fragmentFiles finished.');


function timeintervals = findTimeIntervals
% Returns an array with a row for each time interval that contains selected
% trials.  The first column is the starting timestamp and the second column
% is the ending timestamp.  The starting timestamp is simply the timestamp
% of the trial's lfp_NominalTrialStart.  We want to include the entire BD
% section plus any intertrial interval that may exist between BDOutputStop
% and the next lfp_NominalTrialStart, but not the next
% lfp_NominalTrialStart itself.  However, we DO want to have some sort of
% indication that the recorded segment ends at that time so we don't have
% to rely on the CSC files to infer it.  Therefore, the end time is the
% next trial's lfp_NominalTrialStart timestamp, and we will insert a new
% event ID in place of the terminating lfp_NominalTrialStart (see
% createOneFragmentFile in this file for all the gory exceptions to this
% simple description).  If everything functions correctly, the end result
% of this process is that the output fragment contains only intact trials,
% and should not produce warnings of any extra trial starts or trial ends.
% CSC data are copied in integral frames, so there should always be CSC
% data at the start and end of each <timeinterval>, possibly extending up
% to a frame before the start or after the end.  Nothing is done here
% explicitly to expand time intervals to frame boundaries; this is handled
% in createOneFragmentFile.
%
% In the special case of Joey's data, where lfp_NominalTrialStart =
% BDOutputStop (aka 97), it is logically impossible to insert a splice
% marker without breaking the BD section, so in this case we simply omit
% the event, and the splice marker is the fact that there is a significant
% time delay between BDOutputStart and BDoutputStop.  So the end time has
% to be one Neuralynx clock tick before the BDOutputStop (which, after
% moving the BD sections to lfp_TrialParams, is the event right AFTER the
% lfp_NominalTrialEnd). Except: on the last trial in the fragment, we DO
% need to include the BDOutputStop!
%
% If the start timestamp is at or before the preceding end timestamp, then
% we lengthen the preceding interval, rather than adding another row.
global lfp_Events
lfp_declareGlobals;
timeintervals = [];
selectedtrials = find(lfp_SelectedTrials);
if isequal(lfp_NominalTrialStart, 97)
    % The special case when trial start is the same as Joey's BDOutputStop,
    % and so we assume that trial end is BDOutputStart:
    for trial = selectedtrials
        starttime = lfp_Events(lfp_TrialIndex(trial, 1), 1);
        % At this point, the BD parameters have been removed from
        % lfp_Events, leaving only BDOutputStart and BDOutputStop.
        % <endtime> is 1 microsec before the first event after
        % BDOutputStart (which would normally be BDOutputStop), unless
        % there is no such event, in which case endtime is exactly the last
        % event in the session.
        endtime = lfp_Events( ...
            min(lfp_TrialIndex(trial, 2) + 1, size(lfp_Events,1)), 1 );
        if trial ~= selectedtrials(end)
            endtime = endtime - 1;
        end
        timeintervals = [timeintervals; starttime endtime ];
        if (size(timeintervals,1) > 1) ...
                && (timeintervals(end, 1) <= timeintervals(end-1, 2))
            timeintervals(end-1, 2) = timeintervals(end, 2);
            timeintervals(end, :) = [];
        end
    end
else
    % The normal case, where starttime is the trial start and endtime is
    % the next trial start (except on the last trial, where endtime is
    % exactly the last event in the session):
    for trial = selectedtrials
        starttime = lfp_Events(lfp_TrialIndex(trial, 1), 1);
        if trial < size(lfp_TrialIndex, 1)
            endtime = lfp_Events(lfp_TrialIndex(trial+1, 1), 1);
        else
            endtime = lfp_Events(end, 1);
        end
        timeintervals = [timeintervals; starttime endtime ];
        if (size(timeintervals,1) > 1) ...
                && (timeintervals(end, 1) <= timeintervals(end-1, 2))
            timeintervals(end-1, 2) = timeintervals(end, 2);
            timeintervals(end, :) = [];
        end
    end
end

function createFragments(rulename, timeselection, cscfiles, evfilename, ...
    origtrialnums, sessionName, fragNEV, destdir)
% <rulename> is a string, <timeselection> an N x 2 array.
lfp_declareGlobals;
selectedTrialsIdx = find(lfp_SelectedTrials);
if isempty(destdir)
    destdir = fullfile(lfp_DataDir, rulename);
else
    [p, sessionID] = fileparts(lfp_DataDir);
    destdir = fullfile(destdir, sessionID, rulename);
end
if isempty(selectedTrialsIdx)
    warning('lfp_FragmentFiles:skipping', ...
        'No trials selected, skipping creation of %s', ...
        destdir );
else
    switch exist(destdir)
        case 0
            mkdir(destdir);
        case 7
            warning('lfp_fragmentFiles:overwriting', ...
                'Overwriting like-named files in %s', destdir );
        otherwise
            error('lfp_fragmentFiles:fileexists', ...
                ['There is a non-directory file named ' destdir] );
    end
    % Create info file:
    if isempty(origtrialnums)
        origtrialnums = 1:length(lfp_SelectedTrials);
    end
    trialnums = origtrialnums(selectedTrialsIdx);
    save(fullfile(destdir, 'lfp_FragmentFiles.mat'), 'sessionName', 'trialnums');
    % Find & frag the events file(s):
    DataFiles = dir(lfp_DataDir);
    for file = DataFiles'
        if ~strcmp(file.name, '.') && ~strcmp(file.name, '..')
            [pathstr,name,ext] = fileparts(file.name);
            if ~ismember(upper(ext), {'.MAT', '.EVTSAV'})
                if file.bytes == 16384
                    % empty file; unceremoniously ignore
                    continue
                end
            end
            if strcmpi(file.name, evfilename) && strcmpi(ext, '.MAT')
                ext = 'MAT.EVT'; % symbol in this file for .MAT events file
            end
            if ( (fragNEV && strcmpi(ext, '.NEV')) ...
                    || strcmpi(file.name, evfilename) )
                createOneFragmentFile(rulename, destdir, ...
                    name, ext, timeselection);
            end
        end
    end
    % Frag the CSC file(s):
    cscfiles = reshape(cscfiles, 1, []);
    for filename = cscfiles
        [pathstr,name,ext] = fileparts(char(filename));
        createOneFragmentFile(rulename, destdir, ...
            name, ext, timeselection);
    end
end


function createOneFragmentFile(rulename, destdir, name, ext, timeselection)
% <destdir> is an existing subdirectory of lfp_DataDir.  
% <name>, <ext> are the parts of a filename in lfp_DataDir.
% <timeselection> is N x 2 array.
% All .MAT files are assumed to be CSCs.  The "impossible" value of <ext>
% MAT.EVT denotes an events file in .MAT format.  .NCS files are CSCs and
% .NEV files are events files.
%
% In CSC files, we move our starting time one frame plus one sample
% earlier than the trial start time, to be sure of including the starting
% frame - unless we already have that frame, in which case we must NOT
% include it.  
%
% In events files, the first event is NOT supposed to be a BDOutputStop.  So
% if our first event IS a BDoutputStop (and it normally is in Joey data),
% we need to change the event ID to the novel value 95 = 0x805F.  Also,
% EXCEPT in the special case of lfp_NominalTrialStart = BDOutputStop (aka
% 97; i.e. Joey's data), for each selected time interval that ends with an
% lfp_NominalTrialStart event (which is all of them except for the very
% last trial in the source file) we need to change the event ID of the
% terminating lfp_NominalTrialStart to novel event ID 65535 = 0xFFFF;
% timestamp remains the same.

global lfp_fragmentFiles_save_events lfp_fragmentFiles_save_params
lfp_declareGlobals;
clearance = (lfp_SamplesPerFrame + 1) * lfp_SamplePeriod;
Timestamps = cell(1,0); % each element is also a row
Samples = cell(0);  % as returned by Nlx2MatCSC (points x frames)
EventStrings = cell(0,1);   % each element is a row
% Undo stupid tricks played with input file extension:
if strcmp(upper(ext), 'MAT.EVT')   % symbol in this file for .MAT events file
    fullfilename = fullfile(lfp_DataDir, [name '.MAT']);
else
    fullfilename = fullfile(lfp_DataDir, [name ext]);
end
% Do stupid tricks with output file extension:
if strcmp(upper(ext), '.EVTSAV')
    matfilename = fullfile(destdir, [ name ext ]);
else
    matfilename = fullfile(destdir, [ name '.MAT' ]);
end

%
% Read the input file
%
if ~lfp_NoWaitbar
    hWaitBar = waitbar(0, ['reading ' rulename ' ' name], ...
        'Name', 'Creating fragments');
end
lfp_getEvtIDs;
for tsrow = 1:size(timeselection,1)
    % Note that when the input file is a Matlab file, it gets unnecessarily
    % re-read on every iteration; the read can be pulled out of the loop
    % specifically for Matlab files if that proves to be a performance
    % problem.  The way it is now the code is simpler and more local.
    if ~lfp_NoWaitbar
        waitbar(tsrow/size(timeselection,1), hWaitBar);
    end
    switch upper(ext)
        case {'.NCS' '.MAT'}    % CSC file
            % expand on the left if needed (see header comments on this func)
            ncstimes(1) = max(0, timeselection(tsrow, 1) - clearance);
            if tsrow > 1 && ncstimes(1) <= Timestamps{tsrow-1}(end)
                ncstimes(1) = timeselection(tsrow, 1);
            end
            ncstimes(2) = timeselection(tsrow, 2);
            if strcmp(upper(ext), '.NCS')
                [Timestamps{tsrow} Samples{tsrow} header] = ...
                    dg_readCSC(fullfilename, 'mode', 4, ncstimes);
                if length(Timestamps{tsrow}) ~= size(Samples{tsrow},2)
                    error('lfp_fragmentFiles:badNlx', ...
                        'The file %s contains different numbers of frames and timestamps', ...
                        src );
                end
                if tsrow == 1
                    SamplesUnits = 'AD';
                    conversionfactor = 1;
                    if ~exist('UnitlessFileRegexp') || ...
                            isempty(regexpi(fullfilename, UnitlessFileRegexp))
                        for k = 1:length(header)
                            if regexp(header{k}, '^\s*-ADBitVolts\s+')
                                ADBitVoltstr = regexprep(header{k}, ...
                                    '^\s*-ADBitVolts\s+', '');
                                ADBitVolts = str2num(ADBitVoltstr);
                                if isempty(ADBitVolts)
                                    warning(lfp_fragmentFiles:badADBitVolts', ...
                                        'Could not convert number from:\n%s', ...
                                        header{k} );
                                else
                                    conversionfactor = ADBitVolts;
                                    Samples{tsrow} = conversionfactor ...
                                        * Samples{tsrow};
                                    SamplesUnits = 'V';
                                end
                            end
                        end
                    end
                else
                    Samples{tsrow} = conversionfactor * Samples{tsrow};
                end
            else  % .MAT file
                load('-mat', fullfilename);
                frames2keep = find(...
                    (dg_Nlx2Mat_Timestamps >= ncstimes(1)) ...
                    & (dg_Nlx2Mat_Timestamps <= ncstimes(2)) );
                Timestamps{tsrow} = reshape(...
                    dg_Nlx2Mat_Timestamps(frames2keep), 1, []);
                clear dg_Nlx2Mat_Timestamps;
                Samples{tsrow} = dg_Nlx2Mat_Samples(:,frames2keep);
                clear dg_Nlx2Mat_Samples;
                if exist('dg_Nlx2Mat_SamplesUnits') == 1
                    SamplesUnits = dg_Nlx2Mat_SamplesUnits;
                    clear dg_Nlx2Mat_SamplesUnits;
                else
                    SamplesUnits = 'arbs';
                end
            end
        case {'.NEV' 'MAT.EVT'} % symbol in this file for .MAT events file
            % raw Events file
            if strcmp(upper(ext), '.NEV')
                try
                    [Timestamps{tsrow}, TTL, EventStrings{tsrow,1}] = ...
                        dg_readEvents(fullfilename, 'mode', 4, ...
                        timeselection(tsrow, :) ); %#ok<ASGLU>
                catch e
                    if exist(fullfilename, 'file')
                        Timestamps{tsrow} = [];
                        TTL = [];
                        EventStrings{tsrow,1} = '';
                    else
                        error('lfp_fragmentFiles:nosuchfile', ...
                            'The file %s does not exist', ...
                            fullfilename);
                    end
                end
            else  % .MAT file
                load('-mat', fullfilename);
                evts2keep = find(...
                    (dg_Nlx2Mat_Timestamps >= timeselection(tsrow, 1)) ...
                    & (dg_Nlx2Mat_Timestamps <= timeselection(tsrow, 2)) );
                Timestamps{tsrow} = reshape(...
                    dg_Nlx2Mat_Timestamps(evts2keep), 1, [] );
                clear dg_Nlx2Mat_Timestamps;
                EventStrings{tsrow,1} = dg_Nlx2Mat_EventStrings(evts2keep,:);
                clear dg_Nlx2Mat_EventStrings;
            end
            if tsrow<=size(EventStrings,1) && ~isempty(EventStrings{tsrow,1})
                EventStrings{tsrow,1} = char(EventStrings{tsrow,1});
            end
            if  ~ismember(97, lfp_NominalTrialStart) ...
                    && tsrow<=size(EventStrings,1) ...
                    && ~isempty(EventStrings{tsrow,1})
                % In the normal case, where the last included event is a
                % Trial Start, we substitute ID 0xFFFF for Trial Start: 
                if length(EventStrings) > 2 && isequal(upper(EventStrings{tsrow,1}(end,end-3:end)), ...
                        sprintf('04X', lfp_NominalTrialStart) )
                    EventStrings{tsrow,1}(end,:) = ...
                        'RecID: 4098 Port: 0 TTL Value: 0xFFFF';
                end
            end
        case '.EVTSAV'
            if isempty(lfp_fragmentFiles_save_events)
                % Fully-processed events and trial params
                load('-mat', fullfilename);
                lfp_fragmentFiles_save_events = ...
                    [round(1e6*lfp_save_events(:,1)) lfp_save_events(:,2)];
                lfp_fragmentFiles_save_params = ...
                    lfp_save_params;
            end
            clear lfp_save_events;
            clear lfp_save_params;
            evts2keep = find(...
                (lfp_fragmentFiles_save_events(:,1) >= timeselection(tsrow, 1)) ...
                & (lfp_fragmentFiles_save_events(:,1) <= timeselection(tsrow, 2)) );
            Timestamps{tsrow} = reshape(...
                lfp_fragmentFiles_save_events(evts2keep,1), 1, [] );
            % In this case, EventStrings is actually numeric:
            EventStrings{tsrow,1} = lfp_fragmentFiles_save_events(evts2keep,2);
            % Substitute new ID for terminating Trial Start:
            if isequal(lfp_NominalTrialStart, 97)  && ...
                    EventStrings{tsrow,1}(end) == lfp_NominalTrialStart
                EventStrings{tsrow,1}(end) = 65535;
            end
    end
end

%
% Write the output file
%
if ~lfp_NoWaitbar
    waitbar(tsrow/size(timeselection,1), hWaitBar, ...
        ['writing ' rulename ' ' name]);
end
dg_Nlx2Mat_Timestamps = cell2mat(Timestamps);
if isempty(dg_Nlx2Mat_Timestamps)
    warning('lfp_fragmentFiles:nodata', ...
        'The rule %s did not match any data in file "%s"', ...
        rulename, [name ext] );
end
switch upper(ext)
    case {'.NCS' '.MAT'}
        % We rely on Samples being a vector of frame-formatted
        % matrices:
        dg_Nlx2Mat_Samples = cell2mat(reshape(Samples, 1, []));
        dg_Nlx2Mat_SamplesUnits = SamplesUnits;
        save(matfilename, 'dg_Nlx2Mat_Timestamps', ...
            'dg_Nlx2Mat_Samples', 'dg_Nlx2Mat_SamplesUnits', '-v7.3');
    case {'.NEV' 'MAT.EVT'} % symbol in this file for .MAT events file
        dg_Nlx2Mat_EventStrings = cell2mat(EventStrings);
        % Substitute "dummy" event for beginning BDoutputStop:
        if ~isempty(dg_Nlx2Mat_EventStrings) && ...
                strcmp(dg_Nlx2Mat_EventStrings(1,:), ...
                'RecID: 4098 Port: 0 TTL Value: 0x8061' )
            dg_Nlx2Mat_EventStrings(1,:) = ...
                'RecID: 4098 Port: 0 TTL Value: 0x805F';
        end
        toks = regexp(cellstr(dg_Nlx2Mat_EventStrings), ...
            '0x([\dA-F][\dA-F][\dA-F][\dA-F])', 'tokens');
        TTLstrings = cell(size(toks));
        for toknum = 1:length(toks)
            TTLstrings(toknum) = toks{toknum}{1};
        end
        dg_Nlx2Mat_TTL = hex2dec(TTLstrings);
        save(matfilename, 'dg_Nlx2Mat_Timestamps', ...
            'dg_Nlx2Mat_TTL', 'dg_Nlx2Mat_EventStrings');
    case '.EVTSAV'
        % In this case, EventStrings is actually numeric:
        lfp_save_events = ...
            [ 1e-6*dg_Nlx2Mat_Timestamps' cell2mat(EventStrings) ];
        if ~isempty(lfp_fragmentFiles_save_params)
            lfp_save_params = lfp_fragmentFiles_save_params(find(lfp_SelectedTrials));
        else
            lfp_save_params = lfp_fragmentFiles_save_params;
        end
        save(matfilename, 'lfp_save_events', 'lfp_save_params');
end
if ~lfp_NoWaitbar
    close(hWaitBar);
end

