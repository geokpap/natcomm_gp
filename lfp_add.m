function lfp_add(varargin)
%LFP_ADD calls lfp_fileselect GUI to read additional files after lfp_read
%  has finished executing.
%lfp_add
%lfp_add(..., 'cal')
%lfp_add(..., 'preset', presetDir, presetFiles, presetFormat, newsession)
%lfp_add(..., 'noOffset')
%lfp_add(..., 'sameSession')
%lfp_add(..., 'silent')
%lfp_add(..., 'samperiod')

%OPTIONS
%   lfp_add(..., 'cal')
%       For 'Neuralynx Video (*.NVT, *.DAT)' data, in addition to applying
%       the usual mask for deleting impossible points, a second mask is
%       applied to delete both points from any successive pair that differ
%       by more than one pixel in either dimension.  The remaining points
%       are thus confined to the stationary calibration positions.
%   lfp_add(..., 'deglitch', deglitch_func)
%       <deglitch_func> is a function handle to a function that accepts two
%       input arguments and returns two output arguments which is used for
%       deglitching Neuralynx video tracket data when adding a 'Neuralynx
%       Video (*.NVT, *.DAT)' file.  On both input and output, the first
%       argument is an array of x coordinates, the second is y.
%   lfp_add(..., 'filecheck', funch)
%       Invokes the function handle <funch> in place of the default
%       string-matching procedure to verify the correct placement of
%       channels in an added session in lfp_Samples and lfp_FileNames.
%       <funch> must be a function that accepts a cell string array <X> as
%       its first arg and a string <Y> as its second arg, and returns
%       <true> for each element of <X> that "matches" <Y>.  The default
%       value of <funch> is <@ismember>.  Raises error 'lfp_add:getfilenum'
%       if there is more than one match.
%   lfp_add(..., 'maskParameters',maskParameters)
%       This overrides the default video mask parameters. maskParameters is
%       a struct that has the following fields (see 'maskParameters' in
%       code for default values):
%           maskLabel
%           X0
%           X2
%           X3
%           Y0
%           Y1
%           Y4
%           Y3
%       maskParameters is used as a mask to remove bogus video tracker
%       points from raw data.  The symbols are the same as in "T Maze
%       Layout Parameters.doc", but the values must cover a considerably
%       larger area (i.e. all locations where the rat's head LED could
%       legitimately be).  11-Nov-2009 DG updated the defaults to work with
%       the new (03-Nov-2009) Cheetah computer and video box.  The old
%       values are in lfp_lib in the file 'old_data_default_vt_mask.mat'.
%       To use them, do:
%           >> load('old_data_default_vt_mask.mat')
%           >> lfp_add('maskParameters',maskParameters)
%       Alternatively, for greater generality, maskParameters can have just
%       two fields, namely 'maskLabel' and 'rects', where
%       maskParameters.rects is a four-column array containing two x
%       coordinates and two y coordinates that demark one rectangle on each
%       row.  Col 1 is the left-hand x boundary, col 2 right-hand x
%       boundary, col 3 numerically lower y boundary (visually upper
%       boundary in image coordinates), col 4 numerically upper y boundary
%       (visually lower boundary in image coordinates).  Any data points
%       that are on or within the boundary on ANY row are masked out.  In
%       other words, the masked-out region is the OR of all the rectangles.
%       The values Inf and -Inf can be used and are interpreted as you
%       might hope.  The value [0 0 0 0] can be used to mask out only those
%       usually-bogus points that are at the origin.
%   lfp_add(..., 'noOffset')
%       Works only when adding new sessions.  Adds the new session with
%       zero timestamp offset.  WARNING: This should only be used to rejoin
%       subfragments that were created by splitting a single fragment (or
%       session) into first and second halves.  The first half must be
%       loaded first, and then the second half must then be added using
%       'noOffset'.
%   lfp_add(..., 'preset', presetDir, presetFiles, presetFormat, newsession)
%       This works as for lfp_read2, i.e. the GUI is bypassed, <presetDir>
%       specifies the data directory, <presetFiles> is a cell array that
%       specifies the files to add (filename w/ ext, no directory),
%       <presetFormat> specifies the file format (which must be an exact
%       match to one of the possible values of
%       lfp_fileselect_FileTypesPopup_string; see "switch
%       char(lfp_fileselect_FileTypesPopup_string)" in code). <newsession>
%       must be true or false, and is the answer to the dialog box 'Do you
%       want to add a new session?'.  If true, then the list of files
%       should start with an events file, as for lfp_read2.  Adding a
%       "new session" is the ONLY way to add CSC files and/or Events files,
%       even if the so-called "new session" is actually the 'sameSession'.
%   lfp_add(..., 'sameSession')
%       Skips the error checking that makes sure the session is not already
%       loaded.  Can be used to load multiple fragments from the same
%       session.  Use 'noOffset' in addition to 'sameSession' to preserve
%       the original timestamps.
%   lfp_add(..., 'samperiod', samperiod)
%       Sets the default sample period to use when creating a new CSC
%       channel if there are no other CSC channels loaded.
%   lfp_add(..., 'silent') - opposite of verbose.


%$Rev: 423 $
%$Date: 2023-10-05 16:26:26 -0400 (Thu, 05 Oct 2023) $
%$Author: dgibson $

lfp_declareGlobals;

argnum = 1;
calflag = false;
deglitchfuncflag = false;
filecheck_funch = @ismember;
presetflag = false;
noOffsetflag = false;
sameSessionflag = false;
mazeParametersflag = false;
readNLXVT_opts = {};
silentflag = false;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'cal'
            calflag = true;
        case 'deglitch'
            deglitchfuncflag = true;
            deglitch_funch = varargin{argnum+1};
            argnum = argnum+1;
        case 'filecheck'
            filecheck_funch = varargin{argnum+1};
            argnum = argnum+1;
        case 'preset'
            if length(varargin) < argnum+4
                error('''preset'' option requires 4 additional arguments');
            end
            presetflag = true;
            presetdir = varargin{argnum+1};
            selectedfiles = varargin{argnum+2};
            lfp_fileselect_FileTypesPopup_string = varargin{argnum+3};
            newflag = varargin{argnum+4};
            argnum = argnum + 4;
            if ~isequal(class(selectedfiles), 'cell')
                error('lfp_add:badarg1', ...
                    '<presetFiles> must be a cell array');
            end
            if ~isequal(class(presetdir), 'char')
                error('lfp_add:badarg2', ...
                    '<presetDir> must be a string');
            elseif exist(presetdir) ~= 7
                error('lfp_add:badarg3', ...
                    'Directory "%s" does not exist', presetdir);
            end
            if ~islogical(newflag)
                if isnumeric(newflag)
                    newflag = logical(newflag);
                else
                    error('lfp_add:badnewflag', ...
                        'For option ''preset'', <newsession> must be a logical or numeric value, not %s.', ...
                        class(newflag));
                end
            end
            lfp_DataDir = presetdir;
        case 'noOffset'
            noOffsetflag = true;
        case 'sameSession'
            sameSessionflag = true;
        case 'samperiod'
            readNLXVT_opts{end+1} = varargin{argnum};
            argnum = argnum + 1;
            readNLXVT_opts{end+1} = varargin{argnum};
        case 'maskParameters'
            mazeParametersflag = true;
            maskParameters = varargin{argnum+1};
            argnum = argnum+1;
        case 'silent'
            silentflag = true;
        otherwise
            error('lfp_add:badoption', '%s', ...
                ['The option "' varargin{argnum} '" is not recognized.'] );
    end
    argnum = argnum + 1;
end

if ~presetflag
    switch questdlg('Do you want to add a new session?')
        case 'Cancel'
            return
        case 'Yes'
            newflag = true;
        case 'No'
            newflag = false;
    end
end

if newflag
    if presetflag
        EVfilename = selectedfiles{1};
        if length(selectedfiles) > 1
            CSCfilelist = selectedfiles(2:end);
        else
            lfp_loadSetup;  % need CSCFileRegexp, CSCFileExt
            LFPDataFiles = dir(lfp_DataDir);
            CSCfilelist = {};
            for entrynum = 1:length(LFPDataFiles)
                file = LFPDataFiles(entrynum);
                filenameCSC = file.name;
                pathnameCSC = fullfile(lfp_DataDir, filenameCSC);
                filenameCSC = upper(filenameCSC);
                % This is where we decide which files to open:
                if regexpi(filenameCSC, [CSCFileRegexp '\.' CSCFileExt])
                    CSCfilelist{end+1, 1} = filenameCSC;
                end
            end
        end
    else
        EVfilename = [];
        CSCfilelist = lfp_fileselect(...
            'csc', 'Please select the CSC files:' );
    end
    sessionnum = createNewSession(CSCfilelist, presetflag, EVfilename, ...
        noOffsetflag, sameSessionflag, silentflag, filecheck_funch);
    if sessionnum == 0
        warning('lfp_add:sessionabort', ...
            'Session creation was aborted.');
    end
    return
else
    sessionnum = length(lfp_SessionNames);
end
ts_offset = lfp_SessionOffsets(sessionnum);

if ~presetflag
    selectedfiles = lfp_fileselect;
end
if ~isempty(selectedfiles)
    if ~presetflag
        if exist('lfp_fileselect.mat') == 2
            load('lfp_fileselect.mat');
        else
            error('lfp_add:nofsfile', ...
                'lfp_fileselect.mat does not exist.' );
        end
    end
    switch char(lfp_fileselect_FileTypesPopup_string)
        case 'CSC (old - *.DAT)'
            readCSC(selectedfiles, 'nlx');
        case 'CSC (new - *.NCS)'
            readCSC(selectedfiles, 'nlx');
        case 'CSC (Matlab - *.MAT)'
            readCSC(selectedfiles, 'mat');
        case 'Single Cut Cluster (*.T)'
            readSpikes('T-file', selectedfiles, silentflag);
        case 'Multiple Cut Cluster (*.MAT)'
            readSpikes('multiMAT', selectedfiles, silentflag);
        case 'Naotaka Clusters (*.DWD)'
            readSpikes('naotaka', selectedfiles, silentflag);
        case 'Neuralynx Video (*.NVT, *.DAT)'
            if mazeParametersflag
                readNLXVT_opts(end+1:end+2) = {'mask', maskParameters};
            end
            if deglitchfuncflag
                readNLXVT_opts(end+1:end+2) = {'deglitch', deglitch_funch};
            end
            readNLXVT(calflag, selectedfiles, readNLXVT_opts{:});
        case 'Rodent Clusters (*.Tnn, *.TTn)'
            readSpikes('rodent', selectedfiles, silentflag);
        case 'Rodent Tracker (*.DAT)'
            readRodentTracker(selectedfiles, silentflag);
        case 'Single Electrode (*.NSE)'
            readSpikes('NSE', selectedfiles, silentflag);
        case 'Single Cut Cluster (*.MAT)'
            readSpikes('MAT', selectedfiles, silentflag);
        otherwise
            warning('lfp_add:badfiletype', ...
                'Cannot add files of type %s', ...
                char(lfp_fileselect_FileTypesPopup_string));
    end
end
end


function sessionnum = createNewSession(CSCfilelist, presetflag, EVfname, ...
    noOffsetflag, sameSessionflag, silentflag, filecheck_funch)
% Many (or perhaps all) of the following are altered in this funcion:
global lfp_UseFileSelect lfp_SessionNames lfp_Events lfp_FileNames ...
    lfp_NoWaitbar lfp_DataDir lfp_SamplesUnits lfp_SamplesPerFrame ...
    lfp_SamplePeriod lfp_Samples lfp_TimeStamps lfp_ManualEyeCalib ...
    lfp_SessionOffsets lfp_SessionFirstTrials ...
    lfp_TrialIndex lfp_OrigTrialNums lfp_EventsDataDirs lfp_TrialParams ...
    lfp_SelectedTrials lfp_ActiveFilenums lfp_NewTrialTarget
% used by lfp_loadSetup:
global lfp_SetupName %#ok<NUSED>

% sessionnum = 0 indicates abort
sessionnum = 0;
lfp_loadSetup;
lfp_getEvtIDs; % just in case the setup file didn't call it already
oldUseFileSelect = lfp_UseFileSelect;
if presetflag
    EVFilename1 = EVfname;
    EVFilename2 = [];
    lfp_UseFileSelect = false;
    read_mode = 'preset';
else
    lfp_UseFileSelect = true;
end
if ~exist('usestrobe', 'var')
    usestrobe = true;
end
[newParams, newEvents, ~, badtrials] = ...
    lfp_readEvents(read_mode, EVFilename1, EVFilename2, ...
    false, silentflag, usestrobe);
if ~isempty(badtrials)
    error('lfp_add:badtrials', ...
        'There are bad trials reported by ''lfp_readEvents''.');
end
lfp_UseFileSelect = oldUseFileSelect;
if isempty(newEvents)
    warning('lfp_add:eventsabort', ...
        'Reading of events file was aborted.');
    return
end
newEvents(:,1) = 1.0e-6 * newEvents(:,1);
newIndex = lfp_createTrialIndex(newEvents, numel(newParams));
if ~isempty(lfp_NewTrialTarget)
    newIndex = lfp_handleNewTrialTarget(newIndex, newEvents, ...
        length(newParams));
end

[origtrialnums, sessionname] = lfp_getSessionInfo;
if ~(noOffsetflag || sameSessionflag) ...
        && ismember(sessionname, lfp_SessionNames)
    warning('lfp_add:duplicatesession', ...
        'There is already a session named "%s" loaded.', ...
        sessionname );
end
sessionnum = length(lfp_SessionNames) + 1;

maxtimestamp = findmaxTS;

% Find nice round number offset that allows at least 1 sec clearance
if noOffsetflag
    ts_offset = 0;
else
    ts_offset = ceil((maxtimestamp + 1) / 1000) * 1000;
end

newEvents(:,1) = newEvents(:,1) + ts_offset;
numoldevents = size(lfp_Events,1);

if isequal(CSCfilelist, {''})
    CSCfilelist = {};
else
    for filenameCSC = reshape(CSCfilelist, 1, [])
        [pathstr,name,ext] = fileparts(char(filenameCSC));
        gotfile = ~cell2mat(dg_mapfunc(@isempty, lfp_FileNames));
        if isempty(getfilenum(filenameCSC, filecheck_funch))
            error('lfp_add:filemismatch', ...
                'There is no CSC file that matches "%s" in previously loaded session', ...
                name );
        end
    end
end
% Now actually read and append the CSC files:
if ~lfp_NoWaitbar
    hWaitBar = waitbar(0, '', 'Name', 'Reading files');
end
timestamps = cell(1, length(lfp_FileNames));
for CSCindex = 1:length(CSCfilelist)
    filenum = getfilenum(CSCfilelist{CSCindex}, filecheck_funch);
    % Read the file into memory:
    pathnameCSC = fullfile(lfp_DataDir, CSCfilelist{CSCindex});
    switch lower(ext)
        case '.mat'
            load(pathnameCSC, '-mat');
            newtimestamps{filenum} = dg_Nlx2Mat_Timestamps;
            clear dg_Nlx2Mat_Timestamps;
            newsamples{filenum} = dg_Nlx2Mat_Samples;
            clear dg_Nlx2Mat_Samples;
            if exist('dg_Nlx2Mat_SamplesUnits') == 1
                SamplesUnits = dg_Nlx2Mat_SamplesUnits;
                clear dg_Nlx2Mat_SamplesUnits;
            else
                SamplesUnits = 'arbs';
            end
        otherwise
            [newtimestamps{filenum}, newsamples{filenum}, header] ...
                = Nlx2MatCSC(pathnameCSC, [1, 0, 0, 0, 1], 1, 1);
            if length(newtimestamps{filenum}) ~= size(newsamples{filenum},2)
                error('lfp_add:badNlx', ...
                    'The file %s contains different numbers of frames and timestamps', ...
                    src );
            end
            SamplesUnits = 'AD';
            conversionfactor = 1;
            if ~exist('UnitlessFileRegexp') || ...
                    isempty(regexpi(pathnameCSC, UnitlessFileRegexp))
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
                            newsamples{filenum} = conversionfactor ...
                                * newsamples{filenum};
                            SamplesUnits = 'V';
                        end
                    end
                end
            end
    end
    if ~isequal(SamplesUnits, lfp_SamplesUnits{filenum})
        error('lfp_add:SamplesUnits', ...
            'Filenum %d (%s) has units %s in new session, but %s in previously loaded session.', ...
            filenum, lfp_FileNames{filenum}, ...
            SamplesUnits, lfp_SamplesUnits{filenum} );
    end
    newtimestamps{filenum} = newtimestamps{filenum} * 1e-6;
    newsampleperiod = median(diff(newtimestamps{filenum})) ...
        / lfp_SamplesPerFrame;
    % maxSamplePeriodMismatch must be smaller than lfp_SamplePeriod/2, or
    % else it will screw up lfp_findRecSegments.  However, for most
    % analytic purposes, it should be a LOT smaller than than anyway:
    maxSamplePeriodMismatch = lfp_SamplePeriod/50;
    if abs(newsampleperiod - lfp_SamplePeriod) > maxSamplePeriodMismatch
        if ~lfp_NoWaitbar
            close(hWaitBar);
        end
        error('lfp_add:sampfreqmismatch', ...
            'Cannot add this session due to sample frequency mismatch.' );
    end
    lfp_log(sprintf('Read CSC filenum %d "%s"', ...
        filenum, fullfile(lfp_DataDir, CSCfilelist{CSCindex}) ));
    if ~lfp_NoWaitbar
        waitbar(CSCindex/length(CSCfilelist), hWaitBar);
    end
end

% Reconcile new timestamps, concatenate onto old ones.
if isempty(newtimestamps)
    error('lfp_add:nodata', ...
        'No data to add');
else
    reconciledNewTS = dg_reconcileFrames(newtimestamps);
end
if isempty(reconciledNewTS)
    error('lfp_add:frames2', ...
        'Frame timestamps are irreconcilable!');
end
lfp_TimeStamps = [ reshape(lfp_TimeStamps, 1, []) ...
    reshape(reconciledNewTS, 1, []) + ts_offset ];

% now concatenate the new sampleses onto the old ones.
for CSCindex = 1:length(CSCfilelist)
    filenum = getfilenum(CSCfilelist{CSCindex}, filecheck_funch);
    lfp_Samples{filenum} = [ lfp_Samples{filenum} ...
        newsamples{filenum}( :, ...
        ismember(newtimestamps{filenum}, reconciledNewTS) ) ];
end

neweyecalibs = lfp_readEyeCalib;
if ~isempty(neweyecalibs)
    % Propagate first new calib back to start of new session:
    lfp_ManualEyeCalib = [ lfp_ManualEyeCalib
        {ts_offset}, neweyecalibs(1,2:end) ];
    % Append new calibs with session offset:
    lfp_ManualEyeCalib = [ lfp_ManualEyeCalib
        num2cell(cell2mat(neweyecalibs(:,1)) + ts_offset), ...
        neweyecalibs(:,2:end) ];
end
if ~lfp_NoWaitbar
    close(hWaitBar);
end

% Update the rest of the globals:
lfp_SessionOffsets(sessionnum) = ts_offset;
lfp_SessionNames{sessionnum} = sessionname;
lfp_SessionFirstTrials(sessionnum) = size(lfp_TrialIndex,1) + 1;
if isempty(origtrialnums)
    origtrialnums = 1:size(newIndex, 1);
end
lfp_OrigTrialNums(end+1:end+length(origtrialnums)) = origtrialnums;
lfp_EventsDataDirs{end+1} = lfp_DataDir;
lfp_Events = [ lfp_Events; newEvents ];
lfp_TrialParams = [ lfp_TrialParams newParams ];
lfp_TrialIndex = [ lfp_TrialIndex(:, 1:2); newIndex + numoldevents ];

if ~isempty(lfp_SamplePeriod)
    lfp_createCSCindices;
end

lfp_SelectedTrials(end+1:end+size(newIndex,1)) = ...
    repmat(true, 1, size(newIndex,1));

% Warn of any filenums that have not been appended to:
shortfilenums = [];
maxnumsamples = 0;
for filenum = lfp_ActiveFilenums
    if numel(lfp_Samples{filenum}) > maxnumsamples
        maxnumsamples = numel(lfp_Samples{filenum});
    end
end
for filenum = lfp_ActiveFilenums
    if numel(lfp_Samples{filenum}) < ...
            maxnumsamples
        shortfilenums = [ shortfilenums filenum ];
    end
end
if ~isempty(shortfilenums)
    message = ['The following files are now shorter than others; ' ...
        'this will cause problems if they are used later:' sprintf('\n') ];
    for shortnum = shortfilenums
        message = [ message sprintf(...
            '%2d. %s\n', shortnum, lfp_FileNames{shortnum} )];
    end
    msgbox(message, 'Short file warning', 'warn');
end

% Declare victory and go home:
if ~silentflag
    fprintf(1, ...
        'Added and indexed session %s\n', ...
        lfp_SessionNames{end} );
end
end


function filenum = getfilenum(filename, filecheck_funch)
global lfp_FileNames
[pathstr,name,ext] = fileparts(filename);
gotfile = ~cell2mat(dg_mapfunc(@isempty, lfp_FileNames));
gotfileNums = find(gotfile);
filenum = gotfileNums( ...
    feval(filecheck_funch, upper(lfp_FileNames(gotfile)), upper(name)) );
if length(filenum) > 1
    error('lfp_add:getfilenum', ...
        'More than one member of <upper(lfp_FileNames)> matches %s.', ...
        upper(name) );
end
end


function readCSC(CSCfilelist, read_mode)
% CSCfilelist is a cell string vector of file names relative to
%   lfp_DataDir.
% read_mode is 'mat' or 'nlx'.
lfp_declareGlobals;
if length(lfp_SessionNames) > 1
    error('lfp_add:sorry', ...
        'No code written yet for adding new channels to multi-sessions.');
end
for CSCfn = 1:length(CSCfilelist)
    filenum = max(lfp_ActiveFilenums) + 1;
    lfp_ActiveFilenums(end+1) = filenum;
    pathnameCSC = fullfile(lfp_DataDir, CSCfilelist{CSCfn});
    switch read_mode
        case 'mat'
            load('-MAT', pathnameCSC);
            timestamps{filenum}= double(dg_Nlx2Mat_Timestamps);
            if ~isequal(timestamps{filenum}, round(lfp_TimeStamps * 1e6))
                lfp_ActiveFilenums(end) = [];
                error('lfp_add:badTS2', ...
                    'Timestamps for %s do not match those already loaded', ...
                    pathnameCSC );
            end
            clear dg_Nlx2Mat_Timestamps;
            lfp_Samples{filenum} = double(dg_Nlx2Mat_Samples);
            clear dg_Nlx2Mat_Samples;
            if exist('dg_Nlx2Mat_SamplesUnits') == 1
                lfp_SamplesUnits{filenum} = dg_Nlx2Mat_SamplesUnits;
            end
        case 'nlx'
            [timestamps{filenum}, lfp_Samples{filenum}, header] ...
                = Nlx2MatCSC(pathnameCSC, [1, 0, 0, 0, 1], 1, 1);
            if length(timestamps{filenum}) ~= ...
                    size(lfp_Samples{filenum},2)
                lfp_ActiveFilenums(end) = [];
                error('lfp_read2:badNlx', ...
                    'The file %s contains different numbers of frames and timestamps', ...
                    src );
            end
            if ~isequal(timestamps{filenum}, round(lfp_TimeStamps * 1e6))
                lfp_ActiveFilenums(end) = [];
                error('lfp_add:badTS3', ...
                    'Timestamps for %s do not match those already loaded', ...
                    pathnameCSC );
            end
            lfp_SamplesUnits{filenum} = 'AD';
            if ~exist('UnitlessFileRegexp') || ...
                    isempty(regexpi(pathnameCSC, UnitlessFileRegexp))
                for k = 1:length(header)
                    if regexp(header{k}, '^\s*-ADBitVolts\s+')
                        ADBitVoltstr = regexprep(header{k}, ...
                            '^\s*-ADBitVolts\s+', '');
                        ADBitVolts = str2num(ADBitVoltstr);
                        if isempty(ADBitVolts)
                            warning('lfp_add:badADBitVolts', ...
                                'Could not convert number from:\n%s', ...
                                header{k} );
                        else
                            lfp_Samples{filenum} = ADBitVolts ...
                                * lfp_Samples{filenum};
                            lfp_SamplesUnits{filenum} = 'V';
                        end
                    end
                end
            end
    end
    [p, name, e] = fileparts(CSCfilelist{CSCfn});
    lfp_FileNames{filenum} = name;
    lfp_SelectedFiles(filenum) = true;
end
end


function readSpikes(filetype, filenames, silentflag)
lfp_declareGlobals;
if ~lfp_NoWaitbar
    hWaitBar = waitbar(0, '', 'Name', 'Reading files');
else
    hWaitBar = []; % dummy value to pass into functions that need it
end
for fileindex = 1:length(filenames)
    % For rodent clusters, each file may represent a different session, so
    % this must be inside the fileindex loop:
    fname = filenames{fileindex};
    [origtrialnums, sessionname] = lfp_getSessionInfo(fname);
    sessionnum = find(ismember(lfp_SessionNames, sessionname));
    if isempty(sessionnum)
        stars = '*******************************************************';
        warningmsg = sprintf( ...
            [ '\n\n\n\n%s\nYou are adding data from session "%s",\n' ...
            'but that session is not loaded.\n%s\n\n\n' ], ...
            stars, sessionname, stars );
        warning('lfp_add:nosuchsession', ...
            '%s', warningmsg );
        % Use the most recently loaded session
        sessionnum = length(lfp_SessionNames);
    elseif length(sessionnum) > 1
        warning('lfp_add:sessionnames', ...
            'There is more than one session named %s loaded.', ...
            sessionname);
        sessionnum = sessionnum(end);
    end
    ts_offset = lfp_SessionOffsets(sessionnum);
    [pathstr,name,ext] = fileparts(fname);
    switch filetype
        case {'NSE' 'T-file' 'MAT'}
            readSingleCluster(ts_offset, filetype, fname, name);
            if ~lfp_NoWaitbar
                waitbar(fileindex/length(filenames), hWaitBar);
            end
        case {'naotaka'}
            readNaotaka(ts_offset, fileindex, length(filenames), ...
                fname, name, ext, hWaitBar);
        case {'rodent'}
            readRodent(ts_offset, fileindex, length(filenames), ...
                fname, name, ext, hWaitBar);
        case { 'multiMAT' }
            if isequal(lfp_SetupType, 'rodent')
                % Plexon OFS output from tetrodes is by default numbered
                % with three-digit numbers starting with 000, so we need to
                % rename consistently with tetrode filenames (rewritten
                % & tested by DG 08-Mazy-2008):
                trodenumidx = regexpi(name, '_\d+$');
                if isempty(trodenumidx)
                    warning('lfp_add:mcname', ...
                        'Unknown rodent multiple cut cluster filename style: %s', ...
                        name);
                else
                    % include the '_' in basename:
                    trodefilebasename = name(1:trodenumidx);
                    trodenum = name(trodenumidx+1:end);
                    if length(trodenum) < 3
                        warning('lfp_add:trodenumber', ...
                            'rodent multiple cut cluster trode # not incremented in "%s"', ...
                            name);
                    else
                        name=sprintf('%s%03.0f', ...
                            trodefilebasename, str2num(trodenum)+1);
                    end
                end
            end
            readMultiCluster(ts_offset, filetype, fname, name, silentflag);
    end
end
if ~lfp_NoWaitbar
    close(hWaitBar);
end
end

function readSingleCluster(ts_offset, filetype, fname, name)
lfp_declareGlobals;
spikechannel = length(lfp_Spikes) + 1;
lfp_SpikeNames{spikechannel} = name;
% Read the file into memory:
switch filetype
    case 'NSE'
        % returns double to s:
        s = 1e-6 * dg_readSpike(fullfile(lfp_DataDir, fname));
    case 'T-file'
        % returns double to s:
        [h, s] = lfp_readTfile(fullfile(lfp_DataDir, fname));
        s = 1e-6 * s;
    case 'MAT'
        % File must contain vector of spike timestamps in seconds as
        % first variable saved (does not care what its name is)
        S = load(fullfile(lfp_DataDir, fname));
        fields = fieldnames(S);
        % round to microseconds, return double to s:
        s = round(1e6 * S.(fields{1})) / 1e6;
        clear S;
    otherwise
        error('Bad filetype to readSingleCluster');
end
lfp_Spikes{spikechannel} = s + ts_offset;
lfp_SpikeDataDirs{spikechannel} = lfp_DataDir;
lfp_log(sprintf('Read spike channel %d "%s"', ...
    spikechannel, fullfile(lfp_DataDir, fname) ));
end


function readMultiCluster(ts_offset, filetype, fname, name, silentflag)
% Also reads single cluster files, but tacks 'C1' onto end of name
lfp_declareGlobals;
% Read the file into memory:
switch filetype
    case 'multiMAT'
        % File should contain array with spike timestamps in seconds in col
        % 1 and cluster ID in col 2 as first variable saved (does not care
        % what its name is)
        S = load(fullfile(lfp_DataDir, fname));
        fields = fieldnames(S);
        % returns double to s:
        s = S.(fields{1});
        if any(fix(s(:,2)) ~= s(:,2))
            if all(fix(s(:,1)) == s(:,1))
                warning('lfp_add:cols', ...
                    'Columns appear to be reversed in %s, switching', ...
                    fullfile(lfp_DataDir, fname) );
                s = s(:, [2 1]);
            else
                error('lfp_add:nonint', ...
                    'Both columns of %s contain non-integers', ...
                    fullfile(lfp_DataDir, fname) );
            end
        end
        clear S;
    otherwise
        error('Bad filetype to readMultiCluster');
end
% round to microseconds
s(:,1) = round(1e6 * s(:,1)) / 1e6;
% s now contains timestamps in col 1 and cluster IDs in col 2
if size(s,2) < 2
    %Just got timestamps, treat as single cluster
    clusters = 1;
    if ~silentflag
        disp('Reading single cluster file');
    end
else
    clusters = unique(s(:,2))';
end
for clustnum = clusters
    if size(s,2) < 2
        clustselect = 1:size(s,1);
    else
        clustselect = find(s(:,2) == clustnum);
    end
    spikechannel = length(lfp_Spikes) + 1;
    lfp_SpikeNames{spikechannel} = [ name 'C' num2str(clustnum) ];
    lfp_Spikes{spikechannel} = s(clustselect,1) + ts_offset;
    lfp_SpikeDataDirs{spikechannel} = lfp_DataDir;
    lfp_log(sprintf('Read spike channel %d "%s"', ...
        spikechannel, fullfile(lfp_DataDir, fname) ));
end
end


function readNaotaka(ts_offset, fileindex, numfiles, fname, name, ext, ...
    hWaitBar)
% If there are exactly 5 values in lfp_TrialParams, then in addition to
% reading spike times this also appends the <condition>, <numSaccades>, and
% <yellow> fields from dj_readDataFile to lfp_TrialParams (in that order),
% so param 6 is condition, 7 is numSaccades, and 8 is yellow. Uses as many
% spike channels as there are clusters in the file.  See doc for
% lfp_SpikeNames in lfp_read2. Unpack the spike data one trial at a time.
% Spike times in *.DWD file are in 1 ms time units relative to
% taskStartTime, which is event ID 1.  In case there is disagreement of
% greater than 1 ms between the first saccade time in the *.DWD file and
% that in lfp_Events, mark the trial as bad in lfp_BadTrials.
lfp_declareGlobals;
if ~lfp_NoWaitbar
    waitbar( (fileindex - 1)/numfiles, hWaitBar, 'reading from disk' );
end
[trialdata, numtrials, fileheader, filename] = ...
    dj_readDataFile(fullfile(lfp_DataDir, fname));
if length(trialdata) > size(lfp_TrialIndex, 1)
    warning('lfp_add:extratrials', ...
        '%s has more trials than session %s', ...
        fullfile(lfp_DataDir, fname), lfp_SessionNames{end} );
end

% Check if following assumption is true
ncs=cell(size(trialdata));
[ncs{:}] = deal(trialdata.numClusters);
nc = cell2mat(ncs);
if ~all(nc == nc(1))
    warning('lfp_add:nclust', ...
        'number of clusters is not constant over trials' );
end
chan0 = length(lfp_Spikes);
% We assume here that the number of clusters remains constant over trials:
for clustnum = 1:trialdata(1).numClusters
    % create new empty cluster
    spikechannel = chan0 + clustnum;
    lfp_Spikes{spikechannel} = [];
    digits = regexpi(name, '[0-9]+$');
    trodenum = str2num(name(digits:end));
    lfp_SpikeNames{spikechannel} = sprintf('SE%dC%d', ...
        trodenum, clustnum );
    lfp_SpikeDataDirs{spikechannel} = lfp_DataDir;
    lfp_log(sprintf('Read spike channel %d "%s"', ...
        spikechannel, fullfile(lfp_DataDir, fname) ));
end
for trial = 1:numtrials
    if trial > size(lfp_TrialIndex,1)
        warning('lfp_add:nosuchtrial2', ...
            'Trial %d does not exist in loaded data', trial );
    else
        if length(lfp_TrialParams{trial}) == 5
            lfp_TrialParams{trial}(end+1) = trialdata(trial).condition;
            lfp_TrialParams{trial}(end+1) = trialdata(trial).numSaccades;
            lfp_TrialParams{trial}(end+1) = trialdata(trial).yellow;
        end
        eventrange = lfp_TrialIndex(trial,1) : lfp_TrialIndex(trial,2);
        trialevents = lfp_Events(eventrange,:);
        startidx = find(trialevents(:,2) == 1);
        if isempty(startidx)
            error('lfp_add:noStartTime', ...
                'There is no taskStartTime event in trial %d', ...
                trial);
        end
        if length(startidx) > 1
            error('lfp_add:multiStartTime', ...
                'There are multiple taskStartTime events in trial %d', ...
                trial);
        end
        saccidx = find(ismember(trialevents(:,2), 32:2:44));
        if ~isempty(saccidx)
            sacconset = round(1000 * ( ...
                trialevents(saccidx(1), 1) ...
                - trialevents(startidx, 1) ));
            if abs(sacconset - ( ...
                    trialdata(trial).saccadeOnsetTime ...
                    - trialdata(trial).taskStartTime )) > 1
                warning('lfp_add:TSmismatch', ...
                    'Timestamps do not match DWD file for trial %d', ...
                    trial );
                lfp_BadTrials(end+1) = trial;
            end
            % Ideally there should be an else clause, but it seems that
            % dj_readDataFile routinely returns extra garbage saccade
            % onsets, so there's no easy way to compare.
        end
        
        % trialdata(trial).taskStartTime is normally 0, but sometimes
        % not
        internal_offset = trialevents(startidx, 1) ...
            + 1e-3 * trialdata(trial).taskStartTime;
        for clustnum = 1:trialdata(1).numClusters
            s = 1e-3 * trialdata(trial).spikeTimes{clustnum} ...
                + internal_offset + ts_offset;
            lfp_Spikes{chan0 + clustnum}(end+1:end+length(s)) = s;
        end
    end
    if ~lfp_NoWaitbar
        waitbar( (fileindex - 1 + trial/numtrials)/numfiles, hWaitBar, ...
            'unpacking trials' );
    end
end
end


function readRodent(ts_offset, fileindex, numfiles, fname, name, ext, hWaitBar)
lfp_declareGlobals;
[FileHeader, TrialData] = ...
    dg_ReadRodentFormat(fullfile(lfp_DataDir, fname));
if isempty(TrialData)
    error('lfp_add:notrials', ...
        'The T file %s contains no trial data.', ...
        fullfile(lfp_DataDir, fname) );
end
lfp_log(sprintf('Getting trialnums for Rodent Cluster file %s', ...
    fname ));
trialnums = lfp_getRodentClusterTrialnums(TrialData);
% Note that length(trialnums) can be less than FileHeader.TSize but
% non-empty if the Rodent Cluster file contains empty trials at the end of
% the session.  Also, elements of trialnums can be NaN, denoting an
% unmatchable (usually fragmentary) trial in the Rodent Cluster data, in
% which case we simply delete it from TrialData.
if isempty(trialnums)
    if ~lfp_NoWaitbar
        close(hWaitBar);
    end
    [msgstr,msgid] = lastwarn;
    error('lfp_add:badTrials', ...
        'Trial matching failed for %s:\n%s\n%s', ...
        fname, msgid, msgstr );
end
if any(isnan(trialnums))
    badRCtrials = find(isnan(trialnums));
    TrialData(badRCtrials) = [];
    trialnums(badRCtrials) = [];
end
for clustnum = 1:FileHeader.CSize
    spikechannel = length(lfp_Spikes) + 1;
    lfp_Spikes{spikechannel} = [];
    digits = regexpi(ext, '[0-9]+$');
    trodenum = str2num(ext(digits:end));
    lfp_SpikeNames{spikechannel} = sprintf('%s-T%dC%d', name, ...
        trodenum, clustnum );
    % All that we want here is the spike times.  Uses as many spike
    % channels as there are clusters in the file.  See doc for
    % lfp_SpikeNames in lfp_read2. Unpack the spike data one trial at a
    % time.  Spike times are in 0.1 ms time units relative to um, I don't
    % know what, so we use the timestamp of a reference event to
    % reconstruct the absolute timestamps.  The reference event may be
    % different on each trial, because there is no event that is guaranteed
    % to be present on all trials.
    for trialidx = 1:length(trialnums)
        s = TrialData(trialidx).spikes(:,clustnum);
        reftime=[];
        for refevent = [11, 10, 14, 12, 13, 31, 38, 21, 22]
            if TrialData(trialidx).events(refevent) == 0
                continue
            end
            YasuoGateTime = 1e-4 * TrialData(trialidx).events(refevent);
            eventrange = lfp_TrialIndex(trialnums(trialidx),1) : lfp_TrialIndex(trialnums(trialidx),2);
            trialevents = lfp_Events(eventrange,:);
            reftime = trialevents( ...
                find(ismember(trialevents(:,2), refevent)), ...
                1 );
            if ~isempty(reftime)
                break
            end
        end
        if isempty(reftime)
            if ~lfp_NoWaitbar
                close(hWaitBar);
            end
            error('lfp_add:noref', ...
                'There is no possible reference event for trial %d', ...
                trialnums(trialidx));
        end
        if length(reftime) > 1
            warning('lfp_add:multiref', ...
                'There are %d reference events %d in trial %d; using last.', ...
                length(reftime), refevent, trialnums(trialidx) );
            reftime = reftime(end);
        end
        internal_offset = reftime - YasuoGateTime;
        % Only *trailing* zeros represent padding in the spikes array
        nonzeroidx = find(s);
        if isempty(nonzeroidx)
            % All spikes are at time zero; we must assume that SSize is
            % correct.
            padidx = TrialData(trialidx).header.SSize(clustnum) + 1;
        else
            padidx = find(s(nonzeroidx(1):end) == 0);
            if length(padidx) > 1
                padidx = padidx(1);
            end
            padidx = padidx + nonzeroidx(1) - 1;
        end
        if isempty(padidx)
            lastspikeidx = length(s);
        else
            lastspikeidx = padidx - 1;
        end
        s = 1e-4 * s(1:lastspikeidx) + internal_offset + ts_offset;
        if numel(s) ~= TrialData(trialidx).header.SSize(clustnum)
            error('lfp_add:badTfile', ...
                ['T-file %s clustnum %d trial %d ' ...
                'contains %d spikes, but the trial header states %d'], ...
                fullfile(lfp_DataDir, fname), clustnum, trialidx, ...
                numel(s), TrialData(trialidx).header.SSize(clustnum) );
        end
        lfp_Spikes{spikechannel}(end+1:end+length(s)) = s;
    end
    lfp_SpikeDataDirs{spikechannel} = lfp_DataDir;
    lfp_log(sprintf('Read spike channel %d "%s"', ...
        spikechannel, fullfile(lfp_DataDir, fname) ));
    if ~lfp_NoWaitbar
        waitbar( (fileindex - 1 + clustnum/FileHeader.CSize) ...
            / numfiles, hWaitBar );
    end
end
end


function readRodentTracker(filenames, silentflag)
% This is simply added to the most recently loaded session, after checking
% to be sure session names match.
lfp_declareGlobals;
% For mlint's edification:
global lfp_SessionOffsets

hWaitBar = waitbar(0, '', 'Name', 'Converting tracker data to CSC sampled format');
for fileindex = 1:length(filenames)
    fname = char(filenames(fileindex));
    [origtrialnums, sessionname] = lfp_getSessionInfo(fname(2:end));
    if ~strcmp(sessionname, lfp_SessionNames{end})
        close(hWaitBar);
        error('lfp_add:sessionmismatch', ...
            'The rodent tracker file %s is from a different session: %s', ...
            fname, sessionname );
    end
    [pathstr,name,ext] = fileparts(fname);
    % Add filenums if we are on first session, append to existing filenums
    % otherwise:
    if length(lfp_SessionNames) == 1
        if isempty(lfp_ActiveFilenums)
            filenumX = 1;
        else
            filenumX = max(lfp_ActiveFilenums) + 1;
        end
        filenumY = filenumX + 1;
        lfp_ActiveFilenums(end + 1) = filenumX;
        lfp_ActiveFilenums(end + 1) = filenumY;
        lfp_SelectedFiles(end+1:end+2) = true;
        lfp_FileNames{filenumX} = ['trk X'];
        lfp_FileNames{filenumY} = ['trk Y'];
        lfp_Samples{filenumX} = [];
        lfp_Samples{filenumY} = [];
        lfp_SamplesUnits{filenumX} = 'pix';
        lfp_SamplesUnits{filenumY} = 'pix';
    else
        % Multiple sessions are loaded, find existing filenums
        filenumX = find(ismember(lfp_FileNames, ['trk X']));
        filenumY = find(ismember(lfp_FileNames, ['trk Y']));
        if isempty(filenumX) || isempty(filenumY)
            error('lfp_add:nosuchfilenum', ...
                'Could not find filenums for tracker file.' );
        end
    end
    if isequal(lfp_ActiveFilenums, [1 2])
        % There are no CSCs loaded; use default values for:
        lfp_SamplesPerFrame = 512;
        lfp_SamplePeriod = 1e-3;
        startTS = findminTS;
        numframes = ceil( (findmaxTS - startTS) / ...
            (lfp_SamplePeriod * lfp_SamplesPerFrame) );
        lfp_TimeStamps = startTS + (0:numframes-1) * ...
            lfp_SamplePeriod * lfp_SamplesPerFrame;
        lfp_createCSCindices;
    end
    oldlength = numel(lfp_Samples{filenumX});
    oldnumframes = oldlength/lfp_SamplesPerFrame;
    lfp_Samples{filenumX} = [ lfp_Samples{filenumX} ...
        zeros(1, length(lfp_TimeStamps) * lfp_SamplesPerFrame) ];
    lfp_Samples{filenumY} = [ lfp_Samples{filenumY} ...
        zeros(1, length(lfp_TimeStamps) * lfp_SamplesPerFrame) ];
    fid = fopen(fullfile(lfp_DataDir, fname), 'r');
    if fid == -1
        warning('lfp_add:openfailed', ...
            'Could not open rodent tracker file %s', ...
            fullfile(lfp_DataDir, fname) );
    else
        % skip header lines
        line = dg_ReadLn(fid);
        line = dg_ReadLn(fid);
        line = dg_ReadLn(fid);
        % Read the tracker data
        trackerdata = (fscanf(fid, '%i %i %i %i %i', [5 Inf]))';
        trackerTS = trackerdata(:,1) * 1e-4 + lfp_SessionOffsets(end);
        fclose(fid);
    end
    if any(trackerTS(2:end) <= trackerTS(1:end-1))
        warning('lfp_add:badTS2', ...
            'Some timestamps in %s are out of order!', fname);
    end
    if isempty(trackerdata)
        if ~silentflag
            disp(['No data read from ' ...
                fullfile(lfp_DataDir, fname) ]);
        end
    else
        % Fill each CSC frame with interpolated data; ts is timestamp
        % of the frame we are filling, and trkidx (tracker index)
        % points at the first of the two tracker time points between
        % which we are interpolating.
        framenum = 1 + oldnumframes;
        trkidx = dg_binsearch(trackerTS, lfp_TimeStamps(framenum));
        % Pad beginning if needed:
        if trkidx > 1
            % Tracker data starts before CSC data
            trkidx = trkidx - 1;
        else
            % Tracker data starts after CSC data.
            % Initial frame(s) must be filled with padding up to the
            % time point where tracker data start:
            ts = lfp_TimeStamps(framenum);
            while trackerTS(trkidx) > ts
                % Fill beginning part with padding
                timepoints = ts + (0:511)*lfp_SamplePeriod;
                offset = (framenum-1)*lfp_SamplesPerFrame;
                points2fill = find(timepoints <= trackerTS(trkidx));
                lfp_Samples{filenumX}(oldlength+offset+points2fill) = ...
                    trackerdata(1, 2);
                lfp_Samples{filenumY}(oldlength+offset+points2fill) = ...
                    trackerdata(1, 3);
                % Interpolate ending part
                while (trackerTS(trkidx) <= timepoints(end)) ...
                        && (trkidx < length(trackerTS))
                    points2intrp = find(timepoints > trackerTS(trkidx) ...
                        & timepoints <= trackerTS(trkidx+1) );
                    intrpTrkData(filenumX, offset, timepoints, ...
                        points2intrp, 2, trackerdata, trkidx, trackerTS);
                    intrpTrkData(filenumY, offset, timepoints, ...
                        points2intrp, 3, trackerdata, trkidx, trackerTS);
                    trkidx = trkidx + 1;
                end
                if trkidx >= length(trackerTS)
                    % Ran out of trackerdata, must pad rest of frames
                    padEndTrkData(filenumX, trackerTS, trackerdata, 2);
                    padEndTrkData(filenumY, trackerTS, trackerdata, 3);
                    framenum = length(lfp_TimeStamps);
                else
                    framenum = framenum + 1;
                    trkidx = max(1, trkidx - 1);    % Backtrack to start next frame
                    ts = lfp_TimeStamps(framenum);
                end
            end
        end
        % Fill in time to end of CSC data
        while framenum < length(lfp_TimeStamps)
            waitbar((fileindex - 1 + framenum/length(lfp_TimeStamps)) ...
                / length(filenames), hWaitBar );
            ts = lfp_TimeStamps(framenum);
            timepoints = ts + (0:lfp_SamplesPerFrame-1)*lfp_SamplePeriod;
            offset = (framenum-1)*lfp_SamplesPerFrame;
            % interpolate frame # framenum for X and Y channels
            while (trackerTS(trkidx) <= timepoints(end)) ...
                    && (trkidx < length(trackerTS))
                points2intrp = find(timepoints > trackerTS(trkidx) ...
                    & timepoints <= trackerTS(trkidx+1) );
                intrpTrkData(filenumX, offset, timepoints, ...
                    points2intrp, 2, trackerdata, trkidx, trackerTS);
                intrpTrkData(filenumY, offset, timepoints, ...
                    points2intrp, 3, trackerdata, trkidx, trackerTS);
                trkidx = trkidx + 1;
            end
            if trkidx >= length(trackerTS)
                % Ran out of trackerdata, must pad rest of frame(s)
                padEndTrkData(filenumX, trackerTS, trackerdata, 2);
                padEndTrkData(filenumY, trackerTS, trackerdata, 3);
                framenum = length(lfp_TimeStamps);
            else
                framenum = framenum + 1;
                trkidx = trkidx - 1;    % Backtrack to start next frame
            end
        end
    end
    waitbar(fileindex/length(filenames), hWaitBar);
end
close(hWaitBar);
end


function intrpTrkData(filenum, offset, tpoints, points2intrp, ...
    trkcolnum, trkdata, trkidx, trkTS)
% points2intrp can be empty if there is a gap between frames, in which case we
% do nothing.
if isempty(points2intrp) return
end
global lfp_Samples;
m = (trkdata(trkidx+1,trkcolnum)-trkdata(trkidx,trkcolnum)) ...
    / (trkTS(trkidx+1)-trkTS(trkidx)) ;
delta = ( tpoints(points2intrp) ...
    - trkTS(trkidx) ) * m;
lfp_Samples{filenum}(offset + points2intrp) ...
    = trkdata(trkidx,trkcolnum) + delta;
end


function padEndTrkData(filenum, trackerTS, trackerdata, tkcolnum)
lfp_declareGlobals;
lastTrackerPoint = lfp_time2index(trackerTS(end));
lastsample = length(lfp_TimeStamps) * lfp_SamplesPerFrame;
lfp_Samples{filenum}(lastTrackerPoint:lastsample) ...
    = trackerdata(end, tkcolnum);
end



function maxtimestamp = findmaxTS
% Find maximum value in all timestamps
global lfp_SpikeNames lfp_Spikes lfp_SamplesPerFrame lfp_SamplePeriod lfp_Events lfp_TimeStamps
maxtimestamp = 0;
for cluster = 1:length(lfp_SpikeNames)
    maxtimestamp = max(maxtimestamp, lfp_Spikes{cluster}(end));
end
if ~isempty(lfp_TimeStamps)
    maxtimestamp = max(maxtimestamp, ...
        lfp_TimeStamps(end) + lfp_SamplesPerFrame * lfp_SamplePeriod );
end
maxtimestamp = max(maxtimestamp, max(lfp_Events(:,1)));
end


function mintimestamp = findminTS
% Find minimum value in all timestamps
global lfp_SpikeNames lfp_Spikes lfp_SamplesPerFrame lfp_SamplePeriod lfp_Events lfp_TimeStamps
mintimestamp = Inf;
for cluster = 1:length(lfp_SpikeNames)
    mintimestamp = min(mintimestamp, lfp_Spikes{cluster}(end));
end
if ~isempty(lfp_TimeStamps)
    mintimestamp = min(mintimestamp, ...
        lfp_TimeStamps(end) + lfp_SamplesPerFrame * lfp_SamplePeriod );
end
mintimestamp = min(mintimestamp, min(lfp_Events(:,1)));
end


function readNLXVT(calflag,filename,varargin)
%NLXVTREAD - Reads Neuralynx VT files
%readNLXVT(calflag,filename)
%readNLXVT(...'mask', maskParameters)
%readNLXVT(...'deglitch', deglitch_funch)
%readNLXVT(..., 'samperiod', samperiod)

lfp_declareGlobals

deglitch_funch = @lfp_add_default_deglitch_function;
maskParameters = struct('maskLabel','20091008 default mask',...
    'X0',20,'X2',450,'X3',750,'Y0',140,'Y1',340,'Y4',20,'Y3',450);
argnum = 1;
samperiod = 1e-3;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'deglitch'
            argnum = argnum + 1;
            deglitch_funch = varargin{argnum};
        case 'mask'
            argnum = argnum + 1;
            maskParameters = varargin{argnum};
        case 'samperiod'
            argnum = argnum + 1;
            samperiod = varargin{argnum};
            if ~isscalar(samperiod) || ~isnumeric(samperiod)
                error('lfp_add:samperiod', ...
                    '<samperiod> must be a numeric scalar');
            end
        otherwise
            error('lfp_add:badoption', ...
                'The option "%s" is not recognized.', ...
                dg_thing2str(varargin{argnum}) );
    end
    argnum = argnum + 1;
end

% This first part is taken from the testthat.m script.
[ts, x, y] = dg_readVT( fullfile(lfp_DataDir, filename{1}) );
if any(ts(2:end) <= ts(1:end-1))
    warning('lfp_add:badTS', ...
        'Some timestamps in %s are out of order!', filename{1});
end
ts = ts*1e-6;
x = reshape(x, [], 1);
y = reshape(y, [], 1);

% Replace off-maze locations with NaN:
fieldstrs = fieldnames(maskParameters);
if all(ismember( ...
        {'maskLabel', 'X0', 'X2', 'X3', 'Y0', 'Y1', 'Y4', 'Y3'}, ...
        fieldstrs))
    onmaze = (x>maskParameters.X0 & ...
        x<=maskParameters.X2 & ...
        y>maskParameters.Y0 & ...
        y<maskParameters.Y1) | ...
        (x>maskParameters.X2 & ...
        x<maskParameters.X3 & ...
        y>maskParameters.Y4 & ...
        y<maskParameters.Y3);
elseif all(ismember({'maskLabel', 'rects'}, fieldstrs))
    onmaze = true(size(x));
    for rectnum = 1:size(maskParameters.rects, 1)
        onmaze(x>=maskParameters.rects(rectnum,1) & ...
            x<=maskParameters.rects(rectnum,2) & ...
            y>=maskParameters.rects(rectnum,3) & ...
            y<=maskParameters.rects(rectnum,4)) = false;
    end
else
    error('lfp_add:maskParameters', ...
        'The value given for ''maskParameters'' is missing field(s)');
end
x(~onmaze) = NaN;
y(~onmaze) = NaN;

if calflag
    % Remove non-stationary portions of trace from calibration data; do not
    % deglitch or interpolate
    maxpix = 3;
    calmask = abs(diff(x)) > maxpix;
    calmask = abs(diff(y)) > maxpix | calmask;
    calmask = [calmask; false] | [false; calmask];
    x(calmask) = NaN;
    y(calmask) = NaN;
else
    % Deglitch
    [x, y] = deglitch_funch(x, y);
    % Interpolate NaNs pursuant to "the end of the interpolated section
    % must be within <maxdist> city blocks of the beginning, and no longer
    % than <maxpts> points"
    maxpts = 10;
    maxdist = 35;
    startpt = 1;
    while startpt < length(x)
        % <a> and <b> are the good points surrounding the NaNs
        nans = find(isnan(x(startpt:end)));
        if ~isempty(nans)
            a = startpt + nans(1) - 2;
            nonnan = find(~isnan(x( a+1 : min(a+maxpts-1, length(x)) )));
            if a~= 0 && ~isempty(nonnan)
                b = startpt + nans(1) + nonnan(1) - 2;
                if abs(x(b)-x(a)) + abs(y(b)-y(a)) <= maxdist
                    x(a : b) = linspace(x(a), x(b), nonnan(1) + 1);
                    y(a : b) = linspace(y(a), y(b), nonnan(1) + 1);
                end
                startpt = b;
            else
                nonnan = find(~isnan(x(a+1:end)));
                if ~isempty(nonnan)
                    startpt = a + nonnan(1);
                else
                    break;
                end
            end
        else
            break;
        end
    end
end

% Store position traces
if isempty(lfp_ActiveFilenums)
    newXNum = 1;
else
    newXNum = max(lfp_ActiveFilenums) + 1;
end
newYNum = newXNum+1;
lfp_FileNames{newXNum} = 'VTx';
lfp_FileNames{newYNum} = 'VTy';
lfp_ActiveFilenums = unique([lfp_ActiveFilenums newXNum newYNum]);
lfp_SelectedFiles([newXNum newYNum]) = true;

if isequal(lfp_ActiveFilenums, [1 2])
    % There are no CSCs loaded; use default values for:
    lfp_SamplesPerFrame = 512;
    lfp_SamplePeriod = samperiod;
    startTS = findminTS;
    numframes = ceil( (findmaxTS - startTS) / ...
        (lfp_SamplePeriod * lfp_SamplesPerFrame) );
    lfp_TimeStamps = startTS + (0:numframes-1) * ...
        lfp_SamplePeriod * lfp_SamplesPerFrame;
    lfp_createCSCindices;
end

% convert sampling of x, y from VT to CSC format
fullTSarray = repmat(lfp_TimeStamps, lfp_SamplesPerFrame, 1);
fullTSarray = fullTSarray + repmat( ...
    (0:(lfp_SamplesPerFrame-1))'*lfp_SamplePeriod, ...
    1, size(fullTSarray, 2));
lfp_Samples{newXNum} = NaN(size(fullTSarray));
lfp_Samples{newYNum} = NaN(size(fullTSarray));

% Appallingly, interp1 is too dumb to find the points it CAN interpolate
% between, and just returns all NaN if there are ANY NaNs.  So we must
% explicitly find and interpolate each run of non-NaNs:
[nanidx nonnanidx] = dg_findNaNruns(x);
if isempty(nanidx)
    startswithNaN = false;
else
    startswithNaN = (nanidx(1)==1);
end
for k = 1:length(nonnanidx)
    startpt= nonnanidx(k);
    if startswithNaN
        if k < length(nanidx)
            endpt = nanidx(k+1) - 1;
        else
            endpt = length(x);
        end
    else    % x does not start with NaN
        if k <= length(nanidx)
            endpt = nanidx(k) - 1;
        else
            endpt = length(x);
        end
    end
    % There is "no point" in trying to interpolate between a single point:
    if startpt ~= endpt
        % startpt and endpt now point to the terminals of a run of non-NaNs
        % in x and y; we must compute the starting and ending samples of
        % the corresponding segments in lfp_Samples:
        startsample = lfp_time2index(ts(startpt));
        endsample = lfp_time2index(ts(endpt));
        lfp_Samples{newXNum}(startsample:endsample) = ...
            interp1(ts(startpt:endpt), x(startpt:endpt), ...
            reshape(fullTSarray(startsample:endsample), 1, []), ...
            'linear', NaN );
        lfp_Samples{newYNum}(startsample:endsample) = ...
            interp1(ts(startpt:endpt), y(startpt:endpt), ...
            reshape(fullTSarray(startsample:endsample), 1, []), ...
            'linear', NaN );
    end
end
if calflag
    % Insert NaNs to create line breaks at recording breaks:
    lfp_Samples{newXNum}(lfp_RecSegments(:,2)) = NaN;
    lfp_Samples{newYNum}(lfp_RecSegments(:,2)) = NaN;
end
lfp_SamplesUnits{newXNum} = 'pix';
lfp_SamplesUnits{newYNum} = 'pix';
end


function [x, y] = lfp_add_default_deglitch_function(x, y)
% For Classic Rat T Maze tracker data
x2 = dg_superdeglitch(x, 5, 1);
y2 = dg_superdeglitch(y, 5, 1);
x3 = dg_superdeglitch(x2, 5, 5);
y3 = dg_superdeglitch(y2, 5, 5);
x5 = dg_superdeglitch(x3, 15, 5);
y5 = dg_superdeglitch(y3, 15, 5);
x = x5;
y = y5;
end

