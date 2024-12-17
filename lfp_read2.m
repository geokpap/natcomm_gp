function lfp_read2(varargin)
%lfp_read2
%Change date field to make subversion update the properties:
%23-Jul-2012 16:55:26

%$Rev: 418 $
%$Date: 2022-08-06 20:48:17 -0400 (Sat, 06 Aug 2022) $
%$Author: dgibson $

% Read LFP data for subsequent analysis.
%OPTIONS:
%   'alwaysuseVT' - invokes lfp_rodentPreprocess even when reading events
%       from a *.evtsav file. 
%   'delshortframes' - when short frames are encountered, they are deleted
%       and a warning is raised.  Default is to raise an error.
%   'erreventdup' - when reading rodent files, there is a test for
%       duplicate event IDs that could potentially arise from reading
%       non-matching events (events.nev) and VT events (EACQ*) files.  By
%       default, if a duplicate event is found, a warning is raised and
%       processing continues.  The 'erreventdup' option raises an error
%       instead.
%   'fixeventdup', eventlist, events2skip - <eventlist> is a cell array
%       where each entry is a row vector of possible alternative event IDs
%       that should be in the specified order within each trial as
%       lfp_eventAnalysis.  There is no test here for the existence of all
%       events, but if a duplicated event is found, and one of the copies
%       is preceded and followed by events that match the order in
%       <eventlist>, then  the event IDs for all other instances of the
%       event are incremented by 33000 (e.g. extra event 13s becomes event
%       33013s).  <events2skip> is a list of individual event IDs that
%       should be removed from the events list before checking for
%       duplicates.  Any events that normally may occur more than once, or
%       that do not have a fixed sequential position in <eventlist>, should
%       be on this list.
%   'free' - lfp_fakeEvents is called to fake up the whole session as one
%       trial, instead of lfp_readEvents to read the actual events file.
%       When used together with 'preset', it is still necessary to include
%       an event file name as the first element of <files>, but it can just
%       be the empty string as it gets ignored.
%   'noframecheck' - skips code that marks trials bad when they are missing
%       one or more frames of CSC data
%   'preset' - GUIs are bypassed.  Additional params are required to supply
%       the info that would have been given interactively: <datadir>,
%       <files>.  <datadir> is a string containing the absolute
%       pathname of the data directory.  <files> is a cell array whose
%       first element is the name of the events file without the directory,
%       and subsequent elements are names of CSC files to read (also
%       without directories).  If only the events file name is given, then
%       files are read automatically using the CSCFileRegexp, CSCFileExt
%       mechanism.  If the first file name is '', then no CSC files are
%       read.
%   'silent' - opposite of verbose.
% LFP file names must match a regular expression specified in
% 'lfp_setup.m'.  For the rat version, this means filenames of the form
% 'LFP*.DAT', where '*' matches anything.  See lfp_preamble for further
% info on lfp_setup.
% File name matching is not case sensitive.  File numbers are extracted
% from the LFP file names, if possible, as the string of characters in the
% filename starting at the fourth character and ending at the last
% character before '.DAT'.  If that string does not represent a number,
% then the next consecutive number is used (starting at 1 for the first
% file).

% GLOBAL VARIABLES
%   These are listed in ASCII-alphabetic order, i.e. with A-Z preceding
%   a-z.
% lfp_ActiveFilenums: a row vector containing the file numbers that have
%   been read into the system.  File numbers are assigned on a "best
%   effort" basis to make the file numbers match the numerical part (if
%   there is one) at the tail of the filename.  This number-matching
%   mechanism can be disabled simply by setting CSCFileExt to be something
%   that does not match any extension that you actually use.
% lfp_AlignmentRef: the event ID on which to align waveform data for
%   averaging.  May be a scalar or a vector of IDs, in which case the
%   earliest event in the trial that is in the list is used.
% lfp_AlignStyle: : sets style of display of alignment event e.g. in figure
%   titles; valid values are 'name', 'ID'.  Default is 'ID'.  See also
%   lfp_TrialStyle.
% lfp_AutoassignFilenums: if <true>, then an attempt is made to extract a
%   file number from each CSC file name; this is done by finding the
%   filename extension (i.e. ['\.' CSCFileExt]) and extracting the longest
%   possible string of digits that ends immediately before the '.'.
%   If successful, the file is assigned to the filenum extracted from the
%   name.  Note that this mechanism can be defeated by using a filename
%   extension different from that specified in CSCFileExt. If
%   lfp_AutoassignFilenums is <false>, files are assigned to filenums in
%   sequential order as read.  If lfp_AutoassignFilenums is empty, then it
%   is given the default value <true>.
% lfp_BadTrials: a list of trials that are never included in displays or
%   averages regardless of what lfp_SelectedTrials says.  Useful for
%   eliminating trials that are rejected on purely technical grounds, e.g.
%   the signal of interest is completely obliterated by an artifact.  (See
%   functions lfp_findbadtrials and lfp_markbadRecBreaks.)  Initialized by
%   call to 'lfp_readEvents' if <~freeflag>, otherwise initialized empty;
%   thereafter should be modified only by using 'union'.
% lfp_CLimAll: sets pseudocolor scale for all spectrograms.
% lfp_ClickedTrials: used by lfp_disp to accumulate a list of trialnums
%   that have been clicked on in an overlaid ('ovr') display.
%   Automatically reset to [] on each call to lfp_disp(..., 'ovr').  Should
%   not be modified by other functions, but may be read.
% lfp_ClipBoard: a general-purpose persistent global variable that can be
%   used to pass values into otherwise inaccessible places, e.g. in the
%   rule string given to lfp_selectByRule.
% lfp_DataDir: the directory containing the data.
% lfp_DirStructure: a string specifying the directory structure in which
%   the data reside.  'rodent' is the Yasuo structure, where cluster files
%   and video tracker files must have session-based names and reside in the
%   parent directory of the session directory.  The default structure is
%   the type used by monkey people, where directory name includes both
%   session ID and animal ID, and all files are in the one directory.  Note
%   that this does NOT affect the location of rodent auxiliary event files,
%   which are still assumed to be in the parent directory of the session
%   directory.
% lfp_EventColors: a vector of Matlab ColorSpecs that sets the color code
%   for displaying event markers.
% lfp_EventDefaultColor: a Matlab ColorSpec that is used when displaying
%   event markers that have no color specified in the setup.
% lfp_EventDefaultShape: a legitimate value for the 'Marker' property of a
%   Matlab line object, which is used for plotting event markers in raster
%   displays.
% lfp_EventNames: a cell array containing a character string for each event
%   of interest, to use when referring to events in output.  Also used as
%   an indicator of what range of event IDs to expect in the data.
% lfp_Events: a chronological table of event timestamps (first col) and
%   numeric event codes (second col).
% lfp_EventShapes: legitimate values for the 'Marker' property of a
%   Matlab line object, which are used for plotting event markers in raster
%   displays.
% lfp_EventsDataDirs: a list of the directories from which each session's
%   event file was read.
% lfp_ManualEyeCalib: a cell array whose first column contains
%   timestamps in the same units as lfp_TimeStamps, and the second column
%   contains eye calibration structures of the type returned by
%   lfp_getEyeCalib.  The rows are not guaranteed to be sorted.
% lfp_FileDataDirs: has one entry for each filenum, showing the directory
%   from which data were most recently read for that filenum.  The entry is
%   empty for internally-generated filenums.  This is useful only with
%   a single merged session, i.e. where only one session is loaded, and
%   different data files for the same session were read from different
%   directories.
% lfp_HideGUIs: give this the logical value 'true' to automatically
%   suppress display of the GUI window for GUIs that annoyingly still pop
%   up even when using the 'preset' option (e.g. lfp_save) and that accept
%   the option 'visible'. Default is 'false'.
% lfp_FileNames: cell row vector of LFP filenames; may contain empty cells
%   if LFP file numbers are missing in the set of input file names.  The
%   filename extension is not included in lfp_FileNames.
% lfp_FreqLim: sets y-axis scale for plots that have frequency on the y
%   axis.
% lfp_GetSessionNameFromDir: calculate the session name from the actual
%   current directory structure instead of reading it from
%   'lfp_FragmentFiles.mat'.
% lfp_GlobalNames: a cell row vector of strings containing the names of
%   all global variables used in lfp_lib; it is read from 'lfp_globals.mat'
%   by the script 'lfp_declareGlobals.m'.
% lfp_GoodTrialEvent: a vector of event IDs that designate a good trial in
%   a rodent session.  The good trial event is recorded after the trial,
%   separately from the trial data.
% lfp_LogFileName: absolute pathname of the lfp_lib log file.  If not
%   specified in the setup file, then a filename is automatically generated
%   when you run lfp_read, lfp_read2 or lfp_fragmentFiles.  The filename is
%   in the format 'lfp_libNNNNNNNNN.log', where NNNNNNNNN is a 9-digit
%   random number, and if that also does not exist (which is normally the
%   case), the file is created in the current working directory.  If the
%   value of lfp_LogFileName is 'none' (that value is case-sensitive), then
%   no log messages are written.
% lfp_NewTrialEnd: see lfp_NewTrialTarget
% lfp_NewTrialStart: see lfp_NewTrialTarget
% lfp_NewTrialTarget: if not empty, this value identifies a time in each
%   trial at which new trial end and trial start events are to be inserted,
%   using the values of lfp_NewTrialEnd and lfp_NewTrialStart(1) as the
%   event codes.  The insertion is done after parsing out the BD sections
%   using lfp_NominalTrialEnd and lfp_NominalTrialStart to mark trial
%   boundaries, and after constructing an initial lfp_TrialIndex.  The
%   insertion is only done if lfp_Events does not contain any instances of
%   lfp_NewTrialStart. Legitimate values for lfp_NewTrialTarget have the
%   form
%       {<list-of-eventIDs> <offset-in-s>}
%   i.e. a two-element cell array, where the first element is a list of
%   alternative event IDs to use as a time reference, and the second
%   element is the time in seconds relative to the reference at which the
%   lfp_NewTrialEnd event is inserted (lfp_NewTrialStart is inserted one
%   microsecond later).  The first "old trial" runs from the start of the
%   session to the first lfp_NominalTrialEnd; subsequent "old trials" are
%   bounded by two successive lfp_NominalTrialEnd events.  Only the first
%   instance in one "old trial" of any event in <list-of-eventIDs> will
%   trigger the insertion of a lfp_NewTrialEnd and a lfp_NewTrialStart.
%   After the insertion is completed, the values of lfp_NominalTrialStart
%   and lfp_NominalTrialEnd are updated to have the respective values of
%   lfp_NewTrialStart and lfp_NewTrialEnd. If <StopTask> is defined in the
%   setup or evtIDs file and any trial contains one, then the first
%   StopTask in that trial is used as the trial end event. A similar but
%   more elaborate procedure is invoked by defining <taskStopID>.
%   'StopTask' or 'taskStopID' operates only when lfp_NewTrialTarget is in
%   effect; furthermore, 'taskStopID' operates only when new trial
%   starts/stops are actually being inserted.  'taskStopID' is rather
%   specifically tailored to the requirements of Theresa's data, where
%   specific event codes preceding and following 'taskStopID' are used to
%   disambiguate the true 'taskStopID' from the highest code used to
%   represent <searcheventcount>, both of which are 255.  See code for
%   details.
%       NOTE:  this entire mechanism does not work well with "the rodent
%   way of thought" wherein one must locate the good trial markers and then
%   search backwards to find the trial start: the end result includes bad
%   and incomplete trials as the beginnings of the following good trials.
% lfp_NominalTrialEnd: the event ID (integer) that marks the nominal
%   end of a trial.
% lfp_NominalTrialStart: the event (integer) that marks the nominal
%   start of a trial.
% lfp_NoWaitbar: if set true, disables waitbars that are equipped to be
%   disabled.  Persistent.
% lfp_OrigTrialNums: a row vector of trial numbers taken from the original
%   sequential numbering within the session; these may differ from the
%   internal trial numbers used to index e.g. lfp_SelectedTrials if there
%   are multiple sessions or a fragmented session loaded.
% lfp_RecSegments: start and end sample indices of recorded segments; see
%   lfp_TrialRec and functions lfp_createCSCindices and
%   lfp_findRecSegments.
% lfp_SamplePeriod: the time interval between successive CSC points.  If
%   this variable is empty, then there are no CSC data in memory.
% lfp_Samples: cell row vector of Samples arrays returned by Nlx2MatCSC.
%   Each cell contains the Samples read from one file, with one frame per
%   column. Indices match those for lfp_FileNames.
% lfp_SamplesPerFrame: the number of samples in one frame.
% lfp_SamplesUnits: a vector of strings containing the abbreviated name of
%   the units in which the corresponding lfp_Samples values are measured.
%   When reading data from a Neuralynx CSC file, the conversion is
%   explicitly skipped if the absolute pathname matches UnitlessFileRegexp
%   as assigned in the lfp_getEvtIDs_* file.
%   Possible values of lfp_SamplesUnits are:
%       'AD' - raw Analog to Digital Converter output
%       'V' - Volts
%       'arbs' - created by lfp_createWave with unknown units
%       'pix' - pixels
% lfp_SelectedEventIDs: works like lfp_SelectedTrials to disable display of
%   event markers for unselected event IDs.  Is same length as
%   lfp_EventNames by default.
% lfp_SelectedFiles: works like lfp_SelectedTrials to disable analysis of
%   unselected files (channels).
% lfp_SelectedTrials: a logical row vector where each element is true if
%   the corresponding trial number is currently selected.  In most
%   analyses, trials that are not selected are skipped.  Note that this
%   means that if you say e.g. lfp_disp(10) and trial 10 is unselected,
%   then nothing happens at all.  See also lfp_BadTrials.
% lfp_SelectionRule: a string containing the selection rule applied by the
%   last invocation of lfp_selectByRule (which is called by lfp_select).
% lfp_SessionFirstTrials: the internal trial number of the first trial for
%   each loaded session.
% lfp_SessionNames: a string cell vector containing one session name for
%   each session loaded, in load order.
% lfp_SessionOffsets: a vector of doubles containing an offset in seconds
%   for the timestamps of each session loaded; this offset is added to the
%   raw timestamps at load time, and is chosen so that there is no overlap
%   between the timestamps of different sessions.
% lfp_SetupName: to avoid the overhead of opening a file on every
%   invocation of lfp_getEvtIDs, this contains the setup name after the
%   first time the file is read.
% lfp_SetupType: a string that identifies a particular type of setup
%   file; this is used to invoke specialized processing for particular
%   types of data.  Recognized values are:
%           'monkey'
%           'theresa'
%           'naotaka'
%           'rodent'
%           'wheel'
%   However, other values will not generally cause problems, but will
%   also not invoke any protocol-specific processing (e.g. 'simple').
%   Note that this is an independent setting from lfp_SetupName, although
%   in some cases it has the same value as lfp_SetupName.  It is
%   meant to facilitate the process of multiple people (i.e. multiple
%   values of lfp_SetupName) sharing a common configuration (e.g.
%   'rodent') without the need to add each new person's name to the
%   'switch' statements that control pre-processing.
% lfp_SpikeDataDirs: a list of directories from which each channel of spike
%   data were read.
% lfp_SpikeNames: cell row vector of strings identifying the spike data;
%   for T-files this is just the name of the T-file, without the extension;
%   for Rodent files this is of the form <name>-T<n>C<m>, where <name> is
%   the basename of the file, <n> represents electrode number and <m>
%   represents cluster number within that electrode.  In keeping with other
%   Rodent analysis programs, there are no leading zeros in <m> or <n>.
% lfp_Spikes: a cell row vector where each cell contains a vector of
%   spike times; so lfp_Spikes{i}(j) contains the jth timestamp of cluster
%   i.
% lfp_TimeStamps: the canonical set of timestamps, in seconds, shared by
%   all files read.
% lfp_TrialIndex: a table that represents the trial structure of the data.  
%   It contains one row for each trial. The columns are: start event index,
%   end event index, start sample index, end sample index.  
% lfp_TrialParams: a cell row vector, where each cell contains a row
%   vector of trial parameters (e.g. Joey's "Behavioral Data" sections).
%   Does not include delimiting events for data dump, but does include any
%   format codes at the beginning.  Delimiting events are left in
%   lfp_Events.  Time stamps are not included in lfp_TrialParams because
%   they will be practically the same as the timestamps for the delimiting
%   events.
% lfp_TrialRec: a table with one row for each trial, and columns containing
%   start record sample index, end record sample index.
% lfp_TrialStyle: sets style of display of trial selection e.g. in figure
%   titles; valid values are 'rule', 'trialnums'.  See also lfp_AlignStyle.
% lfp_UseFileSelect: if true, file selection is controlled by user via a
%   GUI; otherwise, files are selected automatically by hard-coded logic.
% lfp_XLimAll: sets x-axis scale for plots that have time on the x axis.
% lfp_XTicksOnAll: makes multi-channel plots show X tick marks on all
%   traces if set to true.
% lfp_Xxcov: sets x-axis scale for cross-covariance plots
% lfp_YPwrLim: sets y-axis scale for plots that have spectral power on the
%   y axis.
% lfp_YLimAll: sets y-axis scale for generic wave plots.
% lfp_YLimNum: a cell vector of 1-by-2 matrices that are automatically used
%   to set ylim for the corresponding filenum.  These values override
%   lfp_YLimAll if non-empty.
% lfp_Yxcov: sets y-axis scale for cross-covariance plots


svnrevision = '$Rev: 418 $';
alwaysuseVTflag = false;
delshortframesflag = false;
erreventdupflag = false;
fixeventdup2skip = [];
fixeventduplist = {};
freeflag = false;
noframecheckflag = false;
presetflag = false;
silentflag = false;
argnum = 1;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'alwaysuseVT'
            alwaysuseVTflag = true;
        case 'delshortframes'
            delshortframesflag = true;
        case 'erreventdup'
            erreventdupflag = true;
        case 'fixeventdup'
            argnum = argnum + 1;
            fixeventduplist = varargin{argnum};
            argnum = argnum + 1;
            fixeventdup2skip = varargin{argnum};
        case 'free'
            freeflag = true;
        case 'noframecheck'
            noframecheckflag = true;
        case 'preset'
            presetflag = true;
            argnum = argnum + 1;
            presetDir = varargin{argnum};
            if ~isequal(class(presetDir), 'char')
                error('Using ''preset'', <presetDir> must be a string');
            end
            % make <presetDir> compatible with 'fileparts':
            [p,n] = fileparts(presetDir);
            if isempty(n)
                presetDir = p;
            end
            argnum = argnum + 1;
            presetFiles = varargin{argnum};
            if ~isequal(class(presetFiles), 'cell')
                error('Using ''preset'', <presetFiles> must be a cell array');
            end
        case 'silent'
            silentflag = true;
        otherwise
            error('lfp_read2:badopt', ...
                'Unrecognized option: %s', dg_thing2str(option) );
    end
    argnum = argnum + 1;
end

if ~silentflag
    disp(sprintf(sprintf('Subversion %s', svnrevision(2:end-1))));
end
lfp_preamble;
if isempty(lfp_DataDir)
    lfp_DataDir = pwd;
end
if isempty(lfp_AlignStyle)
    lfp_AlignStyle = 'ID';
end

lfp_log(sprintf('\n\n\t\tSTARTING FROM CLEARED MEMORY\n'));

% Read the event data
if ~freeflag
    if (exist('EVFilename1') ~= 1)
        EVFilename1 = [];
    end
    if (exist('EVFilename2') ~= 1)
        EVFilename2 = [];
    end
    [lfp_TrialParams, lfp_Events, ~, lfp_BadTrials] = ...
        lfp_readEvents(read_mode, EVFilename1, EVFilename2, ...
        alwaysuseVTflag, silentflag, usestrobe);
    if isempty(lfp_Events)
        return
    end
end

readCSCs = true;
if (exist('lfp_useFileSelect') ~= 1 || ~lfp_useFileSelect) ...
        && isequal(newdatadir, 0)
    % This can indicate a total abort or a reading of events file without
    % any CSC file; if the latter, must still do some more processing:
    if ~isempty(lfp_Events)
        readCSCs = false;
    else
        return
    end
end

if readCSCs
    % Read the CSC waveform data

    % construct CSCfilelist (the list of CSC files to read)
    if lfp_UseFileSelect
        % CSCfilelist is already constructed if lfp_fileselect was called.
    else
        if presetflag && ~isempty(presetFiles(2:end))
            if ~isempty(presetFiles{2})
                CSCfilelist = presetFiles(2:end);
            else
                CSCfilelist = {};
            end
        else
            LFPDataFiles = dir(lfp_DataDir);
            CSCfilelist = {};
            for entrynum = 1:length(LFPDataFiles)
                file = LFPDataFiles(entrynum);
                filenameCSC = file.name;
                pathnameCSC = fullfile(lfp_DataDir, filenameCSC);
                if ispc
                    filenameCSC = upper(filenameCSC);
                end
                % This is where we decide which files to open:
                if regexpi(filenameCSC, [CSCFileRegexp CSCFileExt])
                    CSCfilelist{end+1, 1} = filenameCSC;
                end
            end
        end
    end

    if isempty(CSCfilelist)
        if freeflag
            error('lfp_read2:nofreeCSC', ...
                'You must read at least one CSC file when using the ''free'' option');
        end
    else
        % construct list of numbers to give the CSC files
        filenumbers = [];
        CSCfilelist = reshape(CSCfilelist, [], 1);
        for filenameCSC = CSCfilelist'
            filenameCSC = char(filenameCSC);
            if lfp_AutoassignFilenums
                % Extract channel number from filename if possible:
                suffix = regexpi(filenameCSC, ['\.' CSCFileExt]);
                digits = regexpi(filenameCSC, ['[0-9]+\.' CSCFileExt]);
                number = [];
                if ~isempty(digits)
                    numberstring = filenameCSC(digits(1) : suffix(1) - 1);
                    number = str2num(numberstring);
                end
                if isempty(number)
                    filenumbers(end+1) = 0;
                else
                    filenumbers(end+1) = number;
                end
            else
                if isempty(filenumbers)
                    filenumbers = 1;
                else
                    filenumbers(end+1) = filenumbers(end)+1;
                end
            end
        end
        % Fill in legitimate values for the missing filenumbers.  We want to
        % arrange the filenames with all the numbered ones first, followed by
        % non-numbered ones.
        missing = find(filenumbers == 0);
        filenumbers(missing) = max(filenumbers) + 1 : ...
            max(filenumbers) + length(missing);

        % Now actually read the CSC files:
        % Set up wait bar:
        if ~lfp_NoWaitbar
            hWaitBar = waitbar(0, '', 'Name', 'Reading files');
        end
        lfp_ActiveFilenums = [];
        for CSCindex = 1:length(CSCfilelist)
            filenum = filenumbers(CSCindex);
            [pathstr,name,ext] = fileparts(CSCfilelist{CSCindex});
            lfp_FileNames{filenum} = name;
            lfp_ActiveFilenums(end+1) = filenum;
            % Read the file into memory:
            pathnameCSC = fullfile(lfp_DataDir, CSCfilelist{CSCindex});
            if presetflag
                [pathstr,name,ext] = fileparts(pathnameCSC);
                switch lower(ext)
                    case {'.ncs' '.dat'}
                        read_mode = 'nlx';
                    case '.mat'
                        read_mode = 'mat';
                    otherwise
                        error('lfp_read2:badPreset', ...
                            'Preset CSC files must be .NCS or .MAT' );
                end
            end
            switch read_mode
                case 'nlx'
                    [timestamps{filenum}, lfp_Samples{filenum}, header] ...
                        = dg_readCSC(pathnameCSC);
                    if length(timestamps{filenum}) ~= ...
                            size(lfp_Samples{filenum},2)
                        error('lfp_read2:badNlx', ...
                            'The file %s contains different numbers of frames and timestamps', ...
                            src );
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
                                    warning('lfp_read2:badADBitVolts', ...
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
                case 'mat'
                    load('-MAT', pathnameCSC);
                    timestamps{filenum}= double(dg_Nlx2Mat_Timestamps);
                    clear dg_Nlx2Mat_Timestamps;
                    lfp_Samples{filenum} = double(dg_Nlx2Mat_Samples);
                    clear dg_Nlx2Mat_Samples;
                    if exist('dg_Nlx2Mat_SamplesUnits') == 1
                        lfp_SamplesUnits{filenum} = dg_Nlx2Mat_SamplesUnits;
                    end
            end
            lfp_log(sprintf('Read CSC filenum %d "%s"', filenum, pathnameCSC));
            lfp_FileDataDirs{filenum} = lfp_DataDir;
            if ~lfp_NoWaitbar
                waitbar(CSCindex/length(CSCfilelist), hWaitBar);
            end
        end
        if ~lfp_NoWaitbar
            close(hWaitBar);
        end

        lfp_ManualEyeCalib = lfp_readEyeCalib;

        lfp_SelectedFiles(lfp_ActiveFilenums) = true;
        if ~silentflag
            fprintf( ...
                'Successfully opened CSC files for %d channels in directory:\n%s\n', ...
                length(filenumbers), lfp_DataDir );
            for CSCindex = 1:length(CSCfilelist)
                disp(CSCfilelist{CSCindex});
            end
        end
        if length(lfp_YLimNum) < max(lfp_ActiveFilenums)
            lfp_YLimNum{max(lfp_ActiveFilenums)} = [];
        end

        lfp_TimeStamps = lfp_reconcileFrames(timestamps);
        if isempty(lfp_TimeStamps)
            error('lfp_read2:frames2', ...
                'Frame timestamps are irreconcilable! ');
        end

        lfp_SamplesPerFrame = size(lfp_Samples{lfp_ActiveFilenums(1)}, 1);
        %Convert lfp_TimeStamps from usec to seconds:
        lfp_TimeStamps = 1.0e-6 * lfp_TimeStamps;
        % Find the sample period, ignoring timestamp intervals that
        % indicate a break in recording or other trivial defect:
        framedurs = lfp_TimeStamps(2:end)-lfp_TimeStamps(1:end-1);
        minfd = min(framedurs);
        medianfd = median(framedurs);
        % check left tail:
        isshort = framedurs < 0.9 * medianfd;
        if sum(isshort) > 2
            shortframes = find(isshort);
            % short frame at beginning or end does not matter; short frame
            % in middle is big problem:
            if any(~ismember(shortframes, [1 length(framedurs)]))
                if delshortframesflag
                    warning( 'lfp_read2:frames4', ...
                        'Deleting short frames: %s', ...
                        dg_canonicalSeries(shortframes) );
                    for filenum = lfp_ActiveFilenums
                        lfp_Samples{filenum}(:, isshort) = [];
                    end
                    lfp_TimeStamps(isshort) = [];
                else
                    error( 'lfp_read2:frames4', ...
                        'Short frames: %s', ...
                        dg_canonicalSeries(shortframes) );
                end
            end
        end
        % check right tail (up to 100 breaks in recording is reasonable):
        islong = framedurs > 1.1 * medianfd;
        if sum(islong) > 100
            warning('lfp_read2:frames5', ...
                'There are %d breaks in recording.', ...
                sum(islong));
        end
        % Axe the tails and compute median sample period:
        framedurs((framedurs < 0.9 * medianfd) | ...
                (framedurs > 1.1 * medianfd)) = [];
        mediansample = medianfd / lfp_SamplesPerFrame;
        if std(framedurs(framedurs < minfd + 1.01*mediansample)) > 1e-6
            warning('lfp_read2:frames1', ...
                'CSC frames have sub-sample jitter.');
        end
        if std(framedurs) > mediansample || minfd < medianfd - mediansample
            msgstr = 'CSC frames are of excessively variable duration.';
            msgstr = sprintf('%s\nstd(framedurs)/mediansample = %.3g', ...
                msgstr, std(framedurs)/mediansample);
            msgstr = sprintf( '%s\n%d short frames', msgstr, ...
                sum(framedurs < medianfd - mediansample) );
            warning('lfp_read2:frames3', msgstr);
        end
        lfp_SamplePeriod = mediansample;
        
        if ~silentflag
            fprintf(1, 'Sampling frequency = %g Hz\n', 1/lfp_SamplePeriod);
        end
        totalsamples = lfp_SamplesPerFrame * length(lfp_TimeStamps);
        if ~silentflag
            fprintf(1, 'Total number of samples = %g\n', totalsamples);
        end

        % Convert the CSC waveforms to vectors as well:
        for filenum = lfp_ActiveFilenums
            %     lfp_Samples{filenum} = reshape(lfp_Samples{filenum}, 1, []);
            if numel(lfp_Samples{filenum}) ~= totalsamples
                error('lfp_read2:CSCmismatch2', ...
                    ['The data length in ' lfp_FileNames{filenum} ...
                    ' does not match the number of timestamps in ' ...
                    lfp_FileNames{lfp_ActiveFilenums(1)} ]);
            end
        end
    end
end % if readCSCs

if freeflag
    lfp_fakeEvents;
else
    %Convert lfp_Events timestamps from usec to seconds:
    lfp_Events(:,1) = 1.0e-6 * lfp_Events(:,1);
end

% From this point on, all timestamps are in seconds.

if ~isempty(lfp_NewTrialTarget) || ...
        exist('taskStopID', 'var') && any(lfp_Events(:,2)==taskStopID)
    s1 = warning('off', 'lfp_createTrialIndex:extraStarts');
    s2 = warning('off', 'lfp_createTrialIndex:extraEnds');
end
lfp_TrialIndex = lfp_createTrialIndex(lfp_Events);
if ~isempty(lfp_NewTrialTarget) || ...
        exist('taskStopID', 'var') && any(lfp_Events(:,2)==taskStopID)
    warning(s1.state, 'lfp_createTrialIndex:extraStarts');
    warning(s2.state, 'lfp_createTrialIndex:extraEnds');
end
lfp_SelectedTrials = repmat(true, 1, size(lfp_TrialIndex,1));

if isequal(lfp_SetupType, 'rodent')
    % Idiot check to make sure there's no funny business when merging Nlx
    % events file with VT events file
    for trial = 1:size(lfp_TrialIndex,1)
        trialeventIDs = lfp_Events( ...
            lfp_TrialIndex(trial,1) : lfp_TrialIndex(trial,2), 2 );
        if ~isempty(fixeventdup2skip)
            trialeventIDs(ismember(trialeventIDs, fixeventdup2skip)) = [];
        end
        if length(unique(trialeventIDs)) < length(trialeventIDs)
            msg = sprintf( ...
                'There is a duplicated event ID in trial %.0f.', ...
                trial);
            if erreventdupflag
                error('lfp_read2:eventdup', msg);
            else
                warning('lfp_read2:eventdup', msg);
            end
            if ~isempty(fixeventduplist)
                % Attempt to fix all fixable event duplications.
                for k = 2:(length(fixeventduplist)-1)
                    trialevents = lfp_Events(lfp_TrialIndex(trial,1) : ...
                        lfp_TrialIndex(trial,2), :);
                    if ~isempty(fixeventdup2skip)
                        trialevents(ismember(trialevents(:,2), fixeventdup2skip), :) = [];
                    end
                    evtidx = find(ismember(trialevents(:,2), fixeventduplist{k}));
                    if length(evtidx) > 1
                        evtidx(evtidx==1 | evtidx==length(trialevents)) = [];
                        goodorder = ismember(trialevents(evtidx-1, 2), ...
                            fixeventduplist{k-1} ) & ...
                            ismember(trialevents(evtidx+1, 2), ...
                            fixeventduplist{k+1} );
                        if sum(goodorder) == 1
                            lfp_Events(evtidx(~goodorder) ...
                                + lfp_TrialIndex(trial,1) - 1, 2) = ...
                                lfp_Events(evtidx(~goodorder) ...
                                + lfp_TrialIndex(trial,1) - 1, 2) + 33000;
                            warning('lfp_read2:eventdup3', ...
                                'Replaced event IDs in event(s) %s', ...
                                dg_canonicalSeries(evtidx(~goodorder)) );
                        end
                    end
                end
                trialeventIDs = lfp_Events( ...
                    lfp_TrialIndex(trial,1) : lfp_TrialIndex(trial,2), 2 );
                if ~isempty(fixeventdup2skip)
                    trialeventIDs(ismember(trialeventIDs, fixeventdup2skip)) = [];
                end
                if length(unique(trialeventIDs)) < length(trialeventIDs)
                    warning('lfp_read2:eventdup2', ...
                        'Attempt to fix event duplication in trial %d failed.', ...
                        trial );
                end
            end
        end
    end
end

if ~isempty(lfp_NewTrialTarget)
    [lfp_TrialIndex, lfp_Events] = lfp_handleNewTrialTarget(lfp_TrialIndex, ...
        lfp_Events, length(lfp_TrialParams));
    lfp_SelectedTrials = true(1, size(lfp_TrialIndex,1));
end

[origtrialnums, sessionname] = lfp_getSessionInfo;

if isempty(origtrialnums)
    lfp_OrigTrialNums = 1:size(lfp_TrialIndex, 1);
else
    lfp_OrigTrialNums = origtrialnums;
end

lfp_SessionNames{1} = sessionname;
lfp_SessionOffsets(1) = 0;
lfp_SessionFirstTrials(1) = 1;
lfp_EventsDataDirs{1} = lfp_DataDir;
if length(lfp_EventNames) < max(lfp_Events(:, 2))
    lfp_EventNames{max(lfp_Events(:, 2))} = '';
end

if ~isempty(lfp_SamplePeriod)
    lfp_createCSCindices;
    for k = 2:size(lfp_RecSegments,1)
        gapwindow = lfp_index2time( ...
            [ lfp_RecSegments(k-1,2) lfp_RecSegments(k)] );
        gapevts = lfp_Events(:,1) > gapwindow(1) ...
            & lfp_Events(:,1) < gapwindow(2);
        numevts = sum(gapevts);
        if numevts > 0
            warning('lfp_read2:gapevts', ...
                'There are %d events in the "non-recorded" interval %.6f - %.6f s', ...
                numevts, gapwindow(1), gapwindow(2) );
            gapevts2 = lfp_Events(gapevts,1) > ...
                gapwindow(1) + lfp_SamplesPerFrame * lfp_SamplePeriod ...
                & lfp_Events(gapevts,1) < ...
                gapwindow(2) - lfp_SamplesPerFrame * lfp_SamplePeriod;
            if sum(gapevts2) > 0
                warning('lfp_read2:gapevts2', ...
                    'There are %d events that are more than one frame into the "non-recorded" interval %.6f - %.6f s', ...
                    sum(gapevts2), gapwindow(1), gapwindow(2) );
            end
        end
    end
end

% Checks for pathology: trials missing frames
% Note that this *only* looks from trial start to trial end, time outside
% of trial time will *not* be checked for missing frames.
if ~isempty(lfp_ActiveFilenums) && ~noframecheckflag
    tempbadtrials = [];
    for trialnum = 1:size(lfp_TrialIndex,1)
        startframe = fix((lfp_TrialIndex(trialnum,3)-1)/lfp_SamplesPerFrame)+1;
        endframe = fix((lfp_TrialIndex(trialnum,4)-1)/lfp_SamplesPerFrame)+1;
        % way up above near the excessively variable frames warning the median
        % frame duration <medianfd> is calculated.
        framediffs = diff(lfp_TimeStamps(startframe:endframe));
        if ~isempty(find(abs(framediffs-medianfd)>0.05))
            tempbadtrials = [tempbadtrials trialnum];
        end
    end
    if ~isempty(tempbadtrials)
        lfp_BadTrials = union(lfp_BadTrials, tempbadtrials);
        warning( 'lfp_read2:markbad', ...
            '%d trials marked bad due to missing frames',...
            length(tempbadtrials) );
        lfp_log( sprintf('%d trials marked bad due to missing frames: %s', ...
            length(tempbadtrials), num2str(tempbadtrials)) );
    end
end

if presetflag
    lfp_UseFileSelect = old_UseFileSelect;
end

if isempty(lfp_SamplesUnits)
    lfp_SamplesUnits = cell(size(lfp_FileNames));
end

if ~silentflag
    fprintf(1, 'lfp_read finished.\n');
end

% Fix any problems with event properties globals:
maxID = max(lfp_Events(:,2));
eventGlobLengths = [
    length(lfp_EventColors)
    length(lfp_EventNames)
    length(lfp_EventShapes)
    length(lfp_SelectedEventIDs)
    ];
if any(eventGlobLengths) < maxID
    old_EventColors = lfp_EventColors;
    lfp_EventColors = cell(maxID, 1);
    lfp_EventColors(1:length(old_EventColors)) = old_EventColors;
    lfp_EventColors(max(1,(length(old_EventColors)+1)):end) = {''};
    old_EventNames = lfp_EventNames;
    lfp_EventNames = cell(maxID, 1);
    lfp_EventNames(1:length(old_EventNames)) = old_EventNames;
    old_EventShapes = lfp_EventShapes;
    lfp_EventShapes = cell(maxID, 1);
    lfp_EventShapes(1:length(old_EventShapes)) = old_EventShapes;
    lfp_SelectedEventIDs(end+1:length(lfp_EventColors)) = true;
end
