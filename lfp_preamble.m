%script lfp_preamble

% If presetflag, sets newdatadir and lfp_DataDir to presetDir. Otherwise,
% leaves value returned by uigetdir in variable newdatadir; sets
% lfp_DataDir only if newdatadir is of class 'char' (i.e. only if uigetdir
% completed successfully).

% If <usestrobe> is not assigned a value by lfp_loadSetup, then it is given
% the value true.  It is passed into lfp_readEvents.

% Calls lfp_loadSetup, which executes the appropriate lfp_setup_<name>
% script.  lfp_setup_<name> scripts should call lfp_getEvtIDs at some
% point; this is necessary in order to use symbolic constants defined for
% the event IDs in the setup file, but it is also necessary for some older
% code that uses other symbolic constants aside from event IDs, and assumes
% that lfp_loadSetup has the effect of calling lfp_getEvtIDs somewhere
% along the line.  Newer (as of 11-Nov-2007) code that wants non-event
% symbolic constants explicitly calls lfp_getEvtIDs to get them.  This
% distinction was introduced to make it possible to obtain values of
% symbolic constants without altering the global values that are set by the
% setup files.

%$Rev: 360 $
%$Date: 2015-07-27 14:35:08 -0400 (Mon, 27 Jul 2015) $
%$Author: dgibson $

lfp_declareGlobals;
% Save persistent globals if lfp_lib has already been used during this
% Matlab session; if not, then we do NOT want to wipe out the values in the
% persistent globals file, so we do not save.
if ~isempty(lfp_Events)
    lfp_savePersistentGlobals;
end
% This is as close as we can get to clearing the globals, because actually
% clearing them somehow makes them inaccessible from the base workspace:
lfp_initializeGlobals;
lfp_XTicksOnAll = false;

% Initialization (where non-[] values/classes matter)
timestamps = {};
lfp_Samples = {};
lfp_FileNames = {};
lfp_SelectedFiles = true;
lfp_TrialParams = {};
lfp_Spikes = {};
lfp_SpikeNames = {};
lfp_EventNames = {};
lfp_EventColors = {};
lfp_NoWaitbar = false;

% Get persistent global values:
try
    load('lfp_PersistentGlobalValues.mat');
catch
    lfp_DataDir = '';
end
if ~strcmp(class(lfp_DataDir), 'char')
    lfp_DataDir='';
end

lfp_loadSetup;

if isempty(lfp_NoWaitbar)
    lfp_NoWaitbar = false;
end
if isempty(lfp_LogFileName)
    logname = sprintf('lfp_lib%d.log', round(1e9*rand));
    lfp_LogFileName = which(logname);
    if isempty(lfp_LogFileName)
        lfp_LogFileName = fullfile(pwd, logname);
    end
end
if isempty(lfp_SelectedEventIDs)
    lfp_SelectedEventIDs = logical(ones(1, length(lfp_EventNames)));
end
if length(lfp_EventColors) < length(lfp_EventNames)
    lfp_EventColors{length(lfp_EventNames)} = '';
end
if isempty(lfp_EventDefaultColor)
    lfp_EventDefaultColor = 'b';
end
if isempty(lfp_EventDefaultShape)
    lfp_EventDefaultShape = '.';
end
if isempty(lfp_GoodTrialEvent)
    lfp_GoodTrialEvent = 90;
end
if isempty(lfp_DirStructure)
    lfp_DirStructure = lfp_SetupType;
end
if isempty(lfp_HideGUIs)
    lfp_HideGUIs = false;
end
if isempty(lfp_AutoassignFilenums)
    lfp_AutoassignFilenums = true;
end
if ~exist('usestrobe', 'var')
    usestrobe = true;
end

CSCfilelist = {};
if (exist('useFileSelect') == 1) && logical(useFileSelect)
    lfp_UseFileSelect = true;
else
    lfp_UseFileSelect = false;
end
if (exist('presetflag') == 1) && logical(presetflag)
    old_UseFileSelect = lfp_UseFileSelect;
    lfp_UseFileSelect = false;
    read_mode = 'preset';
    EVFilename1 = presetFiles{1};
    EVFilename2 = [];
end
if lfp_UseFileSelect
    selectedfiles = lfp_fileselect('csc', 'Select CSC files or hit "Cancel" to skip CSCs:');
    newdatadir = 0;
    if ~isempty(selectedfiles)
        newdatadir = lfp_DataDir;
        if exist('lfp_fileselect.mat') == 2
            load('lfp_fileselect.mat');
        else
            error('lfp_preamble:nofsfile', ...
                'lfp_fileselect.mat does not exist.' );
        end
        CSCfilelist = selectedfiles;
        switch char(lfp_fileselect_FileTypesPopup_string)
            case 'CSC (new - *.NCS)'
                read_mode = 'nlx';
            case 'CSC (old - *.DAT)'
                read_mode = 'nlx';
            case 'CSC (Matlab - *.MAT)'
                read_mode = 'mat';
        end
    end
elseif strcmp(read_mode, 'preset')
    newdatadir = presetDir;
    lfp_DataDir = newdatadir;
else
    newdatadir = uigetdir(lfp_DataDir, 'Choose data directory');
    if strcmp(class(newdatadir), 'char')
        lfp_DataDir = newdatadir;
        lfp_savePersistentGlobals;
    else
        return
    end
end
if strcmp(read_mode, 'mat')
    CSCFileExt = 'MAT';
end

