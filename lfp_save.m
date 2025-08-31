function varargout = lfp_save(varargin)
%lfp_save
%lfp_save('double')
%lfp_save('int')
%lfp_save('noclick', taskparams)
%lfp_save('preset', name, ftype)
% Returns [].
% Sets lfp_DataDir as a side-effect.
% The CSC files are exactly like the ones saved by lfp_fragmentFiles.  The
% event files are a new filetype, .EVTSAV (a type of .MAT file), containing
% the value of lfp_Events (with timestamps in seconds), plus the value of
% lfp_TrialParams.  Rodent Clusters are Yasuo-format T-files.
% Multiple Cut Clusters .MAT files are compatible with lfp_add; since these
% files contain absolute cluster ID numbers, an attempt is made to read the
% cluster number from each spike channel's name by looking for the regexp
% pattern 'C\d+$', and if that attempt fails then the cluster is assigned
% an aribtrary sequential ID starting at 1001.
%OPTIONS
%   'destdir', destdir - <destdir> is used in place of lfp_DataDir as the
%       location to which to write the file when using 'preset' option.  If
%       <destdir> does not exist, it is created.
%   'double' - the opposite of 'int', i.e. does not round the values; this
%       is the default.
%   'int' - rounds CSC values to integers, potentially making the saved
%       file roughly up to eight times smaller (depending on the values).
%       Matlab is allowed to make its own decisions how to store the data,
%       so if the data stored as doubles in memory are already integer
%       values, specifying 'int' will not change anything.
%   'noclick' - requires an additional argument, <taskparams>.  Bypasses
%       the GUI for specifying task parameters when saving Rodent Clusters.
%       <taskparams> is the vector of Rodent Clusters file header values:
%       [ ProcType RightCS LeftCS NoGoCS ].
%   'preset' - deprecated.  This option is provided only for backward
%       compatibility.  New code should use lfp_save_noGUI instead of
%       lfp_save('preset'...).
%   'selectedtrials' - save only trials that are not on the lfp_BadTrials
%       list, and for which lfp_SelectedTrials is true.  Works only for
%       saving Rodent Cluster files.
%   'single' - converts values to single-precision float before saving.
%   'visible', val - this option is passed through to Matlab's generic
%       gui_mainfcn to disable display of the GUI window.  <val> must be
%       one of the strings 'on' or 'off' (default is 'on').  Also see
%       lfp_HideGUIs in lfp_read2.

% NOTE FOR FUTURE ENHANCEMENTS: the one sequence of argument types that
% cannot be used is <string>, <handle> (handle being anything that, when
% submitted to ishandle, returns true; see gui_mainfcn).


%      LFP_SAVE, by itself, creates a new LFP_SAVE or raises the existing
%      singleton*.
%
%      H = LFP_SAVE returns the handle to a new LFP_SAVE or the handle to
%      the existing singleton*.
%
%      LFP_SAVE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LFP_SAVE.M with the given input arguments.
%
%      LFP_SAVE('Property','Value',...) creates a new LFP_SAVE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before lfp_save_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to lfp_save_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

%$Rev: 328 $
%$Date: 2014-06-26 13:38:44 -0400 (Thu, 26 Jun 2014) $
%$Author: dgibson $

% Edit the above text to modify the response to help lfp_save

% Last Modified by GUIDE v2.5 30-Aug-2004 13:21:15

global lfp_HideGUIs
if lfp_HideGUIs
    varargin = [{'visible', 'off'}, reshape(varargin(:), 1, [])];
end

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @lfp_save_OpeningFcn, ...
                   'gui_OutputFcn',  @lfp_save_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before lfp_save is made visible.
function lfp_save_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to lfp_save (see VARARGIN)

% Choose default command line output for lfp_save
handles.output = [];

if length(varargin) == 0
    set(handles.MessageLabel, 'String', '');
    showfiletypes = 'all';
end
argnum = 1;
handles.destdir = false;
handles.datatype = 'double';
handles.presetflag = false;
handles.selectedtrialsflag = false;
handles.noclickflag = false;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'destdir'
            argnum = argnum + 1;
            handles.destdir = varargin{argnum};
        case 'double'
             handles.datatype = 'double';
        case 'int'
            handles.datatype = 'int';
        case 'single'
            handles.datatype = 'single';
        case 'preset'
            handles.presetflag = true;
            handles.name = varargin{argnum + 1};
            handles.ftype = varargin{argnum + 2};
            argnum = argnum + 2;
        case 'selectedtrials'
            handles.selectedtrialsflag = true;
        case 'noclick'
            if length(varargin) < argnum + 1
                error('"noclick" option requires array input of form [proctype right_stim left_stim no_go_stim]');
            end
            handles.noclickflag=true;
            handles.inputtaskparams = varargin(argnum + 1);
            argnum = argnum + 1;
        case 'visible'
            if length(varargin) < argnum + 1
                error('"visible" option requires string value of ''on'' or ''off''');
            end
            argnum = argnum + 1;
        otherwise
            error('lfp_save:badoption', ...
                ['The option "' varargin{argnum} '" is not recognized.'] );
    end
    argnum = argnum + 1;
end

% Update handles structure
guidata(hObject, handles);

loadGUIState(handles);
updateDataList(handles);

if handles.presetflag
    % see comments re: 'preset' in this function:
    SaveFileButton_Callback(handles.SaveFileButton, eventdata, handles);
else
    % UIWAIT makes lfp_save wait for user response (see UIRESUME)
    % The really annoying thing about this is that it has to be commented-out
    % in order to debug, but then it has to execute in order for the figure not
    % to automatically close itself.
    uiwait(handles.figure1);
end


function updateDataList(handles)
lfp_declareGlobals;
set(handles.DirLabel, 'String', lfp_DataDir);
value = get(handles.FileTypesPopup, 'Value');
strings = get(handles.FileTypesPopup, 'String');
switch strings{value}
    case 'CSC Channel (*.MAT)'
        set(handles.DataList, 'Value', 1);
        set(handles.DataList, 'String', lfp_FileNames(lfp_ActiveFilenums));
    case 'Events (*.EVTSAV)'
        set(handles.DataList, 'String', 'lfp_Events');
    case 'Rodent Clusters (*.Tnn)'
        set(handles.DataList, 'Value', 1);
        set(handles.DataList, 'String', lfp_SpikeNames);
    case 'Multiple Cut Cluster (*.MAT)'
        set(handles.DataList, 'Value', 1);
        set(handles.DataList, 'String', lfp_SpikeNames);
    otherwise
        error('lfp_save:internalError1', 'programming error');
end


function saveGUIState(handles)
if exist('lfp_save.mat') == 2
    pathname = which('lfp_save.mat');
else
    pathname = 'lfp_save.mat';
end
lfp_save_FileTypesPopup_value = get(handles.FileTypesPopup, 'Value');
strings = get(handles.FileTypesPopup, 'String');
lfp_save_FileTypesPopup_string = ...
    strings(lfp_save_FileTypesPopup_value);
save(pathname, 'lfp_save_FileTypesPopup_value', ...
    'lfp_save_FileTypesPopup_string');


function loadGUIState(handles)
if exist('lfp_save.mat') == 2
    load('lfp_save.mat');
    if lfp_save_FileTypesPopup_value <= length(get(handles.FileTypesPopup, 'String'))
        set(handles.FileTypesPopup, ...
            'Value', lfp_save_FileTypesPopup_value );
    else
        set(handles.FileTypesPopup, 'Value', 1);
    end
end


% --- Outputs from this function are returned to the command line.
function varargout = lfp_save_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
if isempty(handles)
    varargout{1} = {};
else
    varargout{1} = handles.output;
    close(hObject);
end


% --- Executes during object creation, after setting all properties.
function DataList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DataList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in DataList.
function DataList_Callback(hObject, eventdata, handles)
% hObject    handle to DataList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns DataList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from DataList


% --- Executes on button press in SaveFileButton.
function SaveFileButton_Callback(hObject, eventdata, handles)
% hObject    handle to SaveFileButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%   For better or worse, the strategy for handling the 'preset' option is
%   to completely ignore the GUI, which may thus exist invisibly or not
%   exist at all, and to duplicate here whatever operations would have been
%   done to set up the GUI with the proper filenames etc.  The only
%   elements of the GUI that need exist are handles.name and handles.ftype,
%   which are used only for processing 'preset'.
lfp_declareGlobals;
if handles.presetflag
    switch lower(handles.ftype)
        case 'evt'
            filetype = 'Events (*.EVTSAV)';
        case 'csc'
            filetype = 'CSC Channel (*.MAT)';
        case 'rod'
            filetype = 'Rodent Clusters (*.Tnn)';
        case 'mcc'
            filetype = 'Multiple Cut Cluster (*.MAT)';
        otherwise
            error('lfp_save:badpresetext', ...
                'The ftype given with ''preset'' option is not recognized' );
    end
else
    typenum = get(handles.FileTypesPopup, 'Value');
    typestrings = get(handles.FileTypesPopup, 'String');
    filetype = typestrings{typenum};
end
if isequal(handles.destdir, false)
    OutputPathName = lfp_DataDir;
else
    OutputPathName = handles.destdir;
    if ~isempty(OutputPathName) && ~exist(OutputPathName, 'dir')
        mkdir(OutputPathName);
    end
end
set(handles.SaveFileButton, 'Visible', 'off');
set(handles.ChDirButton, 'Visible', 'off');
set(handles.CancelButton, 'Visible', 'off');
switch filetype
    case 'CSC Channel (*.MAT)'
        if handles.presetflag
            channelnames = lfp_FileNames(lfp_ActiveFilenums);
            channels = find(ismember(channelnames, handles.name));
            if isempty(channels)
                error('lfp_save:badpresetname', ...
                    'The name given with the preset option must match a loaded CSC channel' );
            end
        else
            channels = get(handles.DataList, 'Value');
        end
        dg_Nlx2Mat_Timestamps = round(1e6*lfp_TimeStamps);
        disp(sprintf('Saving filenums %s', ...
            dg_thing2str(lfp_ActiveFilenums(channels))));
        for filenum = lfp_ActiveFilenums(reshape(channels, 1, []))
            if handles.presetflag
                OutputFileName = [lfp_FileNames{filenum} '.mat'];
            else
                [OutputFileName, OutputPathName] = ...
                    uiputfile( '*.mat', ...
                    ['Saving ' lfp_FileNames{filenum}], ...
                    fullfile(OutputPathName, ...
                    [lfp_FileNames{filenum} '.mat'] ) );
                if isequal(OutputFileName, 0)
                    CancelButton_Callback([], [], handles);
                    return
                end
            end
            set(handles.MessageLabel, 'String', sprintf(...
                'Saving %s', OutputFileName ));
            refresh(handles.figure1);
            extrasamples = rem(numel(lfp_Samples{filenum}), lfp_SamplesPerFrame);
            padding = [];
            if extrasamples ~= 0
                msgbox(sprintf(...
                    'NOTE: padded "%s" with zeros to fill last frame', OutputFileName ));
                padding = zeros(1, lfp_SamplesPerFrame-extrasamples);
            end
            dg_Nlx2Mat_Samples = reshape(...
                [lfp_Samples{filenum}(1:end), padding], ...
                lfp_SamplesPerFrame, [] );
            switch handles.datatype
                case 'int'
                dg_Nlx2Mat_Samples = round(dg_Nlx2Mat_Samples);
                case 'single'
                dg_Nlx2Mat_Samples = single(dg_Nlx2Mat_Samples);
                case 'double'
                % do nothing
                otherwise
                    error('Doh!');
            end
            dg_Nlx2Mat_SamplesUnits = lfp_SamplesUnits{filenum};
            save(fullfile(OutputPathName, OutputFileName), ...
                'dg_Nlx2Mat_Timestamps', 'dg_Nlx2Mat_Samples', ...
                'dg_Nlx2Mat_SamplesUnits', '-v7.3' );
            lfp_log(sprintf('Saved CSC file %s', fullfile(OutputPathName, OutputFileName)));
        end
    case 'Events (*.EVTSAV)'
        if handles.presetflag
            OutputFileName = [handles.name '.EVTSAV'];
        else
            [OutputFileName, OutputPathName] = ...
                uiputfile(fullfile(OutputPathName, ...
                ['Events.EVTSAV'] ), ...
                'Saving Events' );
            if isequal(OutputFileName, 0)
                CancelButton_Callback([], [], handles);
                return
            end
        end
        set(handles.MessageLabel, 'String', sprintf(...
            'Saving %s', OutputFileName ));
        refresh(handles.figure1);
        lfp_save_events = lfp_Events;
        lfp_save_params = lfp_TrialParams;
        save(fullfile(OutputPathName, OutputFileName), ...
            'lfp_save_events', 'lfp_save_params');
        lfp_log(sprintf('Saved Events file %s', fullfile(OutputPathName, OutputFileName)));
    case 'Rodent Clusters (*.Tnn)'
        % Because all time stamps are saved relative to the estimated start
        % of recording for each trial, this "should work" equally well on
        % any session when multiple sessions are loaded.
        if handles.presetflag
            channelnames = lfp_SpikeNames;
            channels = find(ismember(channelnames, handles.name));
            if isempty(channels)
                error('lfp_save:badpresetname2', ...
                    'The name given with the preset option must match a loaded spike channel' );
            end
        else
            channels = get(handles.DataList, 'Value');
        end
        disp(sprintf('Saving spike channels %s', dg_thing2str(channels)));
        OutputFileName = makeClusterFilename('t', ...
            lfp_SpikeNames{channels(1)});
        if ~handles.presetflag
            [OutputFileName, OutputPathName] = ...
                uiputfile(fullfile(OutputPathName, OutputFileName), ...
                'Saving Rodent Clusters' );
            if isequal(OutputFileName, 0)
                CancelButton_Callback([], [], handles);
                return
            end
        end
        set(handles.MessageLabel, 'String', sprintf(...
            'Saving %s', OutputFileName ));
        refresh(handles.figure1);
        if handles.selectedtrialsflag
            trials = lfp_enabledTrials(1:size(lfp_TrialIndex,1));
        else
            trials = 1:size(lfp_TrialIndex,1);
        end
        % Create file header:
        FileHeader.Format = hex2dec('FFF3');
        clockval = clock;
        FileHeader.Year     = clockval(1);
        FileHeader.Month    = clockval(2);
        FileHeader.Day      = clockval(3);
        FileHeader.SSize    = 0;
        channels = reshape(channels, 1, []);
        if handles.noclickflag
            taskparams = cell2mat(handles.inputtaskparams);
        else
            taskparams = lfp_getRodentTaskParams;
        end
        FileHeader.CSize    = length(channels);   % number of clusters in file
        FileHeader.TSize    = length(trials);   % number of trials
        FileHeader.ProcType = taskparams(1);
        FileHeader.RightCS  = taskparams(2);
        FileHeader.LeftCS   = taskparams(3);
        FileHeader.NoGoCS   = taskparams(4);
        % Re-package spikes & event data according to trial structure, and
        % save the file.  Note that trial numbers are strictly sequential,
        % and may thus disagree with the trial numbers that are saved by
        % Yasuo's programs.  Spike and event times are converted to
        % tenths of milliseconds relative to a fictitious starting time
        % that is arbitrarily set to coincide with the start of the
        % "pre-roll" period.  The actual nominal trial start time cannot be
        % used because then the start event would have a zero timestamp and
        % appear in the Yasuo format to be missing.
        TrialData = [];
        preRoll = 2.0012;
        postRoll = 0.5;
        for trialidx = 1:length(trials)
            trial = trials(trialidx);
            spikes = lfp_getTrialSpikes(trial, channels, preRoll, postRoll);
            for channelidx = 1:length(channels)
                FileHeader.SSize = ...
                    FileHeader.SSize + length(spikes{channelidx});
            end
            TrialData(trialidx).header.TNumber = trial;
            % 12 is value of dg_WRF_MaxCluster from dg_WriteRodentFormat:
            TrialData(trialidx).header.SSize = zeros(1, 12);
            for cluster = 1:FileHeader.CSize
                TrialData(trialidx).header.SSize(cluster) = ...
                    numel(spikes{cluster});
            end
            TrialData(trialidx).header.TSSize = ...
                sum(TrialData(trialidx).header.SSize);
            TrialData(trialidx).header.Free = [1 1 1];
            paddedspikes = zeros(max(TrialData(trialidx).header.SSize), ...
                length(spikes));
            reftime = lfp_Events(lfp_TrialIndex(trial,1),1) - preRoll;
            for cluster = 1:FileHeader.CSize
                paddedspikes(1:length(spikes{cluster}), cluster) = ...
                    round(1e4*(spikes{cluster} - reftime));
            end
            TrialData(trialidx).spikes = paddedspikes;
            % Create Yasuo-style events array, assuming there is no more
            % than one of each event ID per trial.  Include the event
            % before lfp_NominalTrialStart if it is ID 1 (Record On).
            maxEvent = 50; % value of dg_WRF_MaxEvent from dg_WriteRodentFormat
            eventsarray = zeros(1, maxEvent);
            startevtidx = lfp_TrialIndex(trial,1);
            if startevtidx > 1 && lfp_Events(startevtidx - 1, 2) == 1
                startevtidx = startevtidx - 1;
            end
            for evtidx = startevtidx : lfp_TrialIndex(trial,2)
                eventsarray(lfp_Events(evtidx,2)) = round(...
                    1e4*(lfp_Events(evtidx,1) - reftime) );
            end
            if length(eventsarray) < maxEvent
                eventsarray(end+1:maxEvent) = 0;
            end
            TrialData(trialidx).events = eventsarray(1:maxEvent);
            % These are correct SType values when there is an event [31
            % 38 21 22], based on \\Matrisome2\RData\A75\acq02\EACQ02.DAT:
            if eventsarray(31) || eventsarray(21)
                TrialData(trialidx).header.SType = 1;
            elseif eventsarray(38) || eventsarray(22)
                TrialData(trialidx).header.SType = 8;
            else
                TrialData(trialidx).header.SType = 0;
            end
            % These RType values are correct regarding direction:
            if eventsarray(15)
                % Right Turn
                if TrialData(trialidx).header.SType == FileHeader.RightCS
                    % correct
                    TrialData(trialidx).header.RType = 1;
                else
                    % incorrect
                    TrialData(trialidx).header.RType = 4;
                end
            elseif eventsarray(16)
                % Left Turn
                if TrialData(trialidx).header.SType == FileHeader.LeftCS
                    % correct
                    TrialData(trialidx).header.RType = 2;
                else
                    % incorrect
                    TrialData(trialidx).header.RType = 5;
                end
            else
                % Score any other response as "incomplete"; note this does
                % not cover the no-go task:
                TrialData(trialidx).header.RType = 0;
            end
        end
        dg_WriteRodentFormat(fullfile(OutputPathName, OutputFileName), ...
            FileHeader, TrialData);
        lfp_log(sprintf('Saved Rodent Clusters file %s', ...
            fullfile(OutputPathName, OutputFileName)));
    case 'Multiple Cut Cluster (*.MAT)'
        % For compatibility with lfp_add, file must contain array with
        % spike timestamps in seconds in col 1 and cluster ID in col 2 as
        % first variable saved (does not care what its name is)
        if handles.presetflag
            channelnames = lfp_SpikeNames;
            channels = find(ismember(channelnames, handles.name));
            if isempty(channels)
                error('lfp_save:badpresetname3', ...
                    'The name given with the preset option must match a loaded spike channel' );
            end
        else
            channels = get(handles.DataList, 'Value');
        end
        disp(sprintf('Saving spike channels %s', dg_thing2str(channels)));
        OutputFileName = makeClusterFilename('m', lfp_SpikeNames{channels(1)});
        if ~handles.presetflag
            [OutputFileName, OutputPathName] = ...
                uiputfile(fullfile(OutputPathName, OutputFileName), ...
                'Saving Rodent Clusters' );
            if isequal(OutputFileName, 0)
                CancelButton_Callback([], [], handles);
                return
            end
        end
        set(handles.MessageLabel, 'String', sprintf(...
            'Saving %s', OutputFileName ));
        refresh(handles.figure1);
        numrows = 0;
        for ch = reshape(channels, 1, [])
            numrows = numrows + numel(lfp_Spikes{ch});
        end
        lfp_save_spikes = zeros(numrows, 2);
        nextrow = 1;
        arbseqnum = 1000;
        for ch = reshape(channels, 1, [])
            idx = regexp(lfp_SpikeNames{ch}, 'C\d+$');
            if isempty(idx)
                arbseqnum = arbseqnum + 1;
                clust = arbseqnum;
            else
                clust = str2num(lfp_SpikeNames{ch}(idx+1:end));
            end
            numspikes = numel(lfp_Spikes{ch});
            lfp_save_spikes(nextrow : nextrow + numspikes - 1, :) = ...
                [ reshape(lfp_Spikes{ch}, [], 1) ...
                repmat(clust, numspikes, 1) ];
            nextrow = nextrow + numspikes;
        end
        save(fullfile(OutputPathName, OutputFileName), 'lfp_save_spikes');
        lfp_log(sprintf('Saved Multiple Cut Clusters file %s', ...
            fullfile(OutputPathName, OutputFileName)));
    otherwise
        error('lfp_save:internalError2', 'programming error');
end
if ~handles.presetflag
    uiresume(handles.figure1);
end


% --- Executes on button press in ChDirButton.
function ChDirButton_Callback(hObject, eventdata, handles)
% hObject    handle to ChDirButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
lfp_declareGlobals;
newdatadir = uigetdir(lfp_DataDir, 'Choose data directory');
if strcmp(class(newdatadir), 'char')
    lfp_DataDir = newdatadir;
    lfp_savePersistentGlobals;
    updateDataList(handles);
end


% --- Executes on button press in CancelButton.
function CancelButton_Callback(hObject, eventdata, handles)
% hObject    handle to CancelButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(handles.figure1);

% --- Executes during object creation, after setting all properties.
function FileTypesPopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FileTypesPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in FileTypesPopup.
function FileTypesPopup_Callback(hObject, eventdata, handles)
% hObject    handle to FileTypesPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
updateDataList(handles);
saveGUIState(handles);



function OutputFileName = makeClusterFilename(ftype, spikename)
% Construct spike file name from spike name; <ftype> = 't' for Rodent
% Cluster "T-file", <ftype> = 'm' for 'Multiple Cut Cluster (*.MAT)' file.
switch ftype
    case 'm'
        fmt1 = '%s-T%02d.mat';
        fmt2 = 'SessionID-Tnn.mat';
    case 't'
        fmt1 = '%s.T%02d';
        fmt2 = 'SessionID.Tnn';
    otherwise
        error('makeClusterFilename:badtype', 'Bad ftype: %s', ftype);
end
dashindex = strfind(spikename, '-');
Tindex = strfind(spikename, 'T');
Cindex = strfind(spikename, 'C');
underscoreindex = strfind(spikename, '_');
if underscoreindex
    sessionID = spikename(1 : underscoreindex - 1);
    trodenum = str2num( ...
        spikename(underscoreindex + 1 : Cindex - 1) ) ;
    OutputFileName = sprintf(fmt1, ...
        sessionID, ...
        trodenum );
elseif isempty(dashindex) || isempty(Tindex) || isempty(Cindex)
    OutputFileName = fmt2;
else
    trodenum = [];
    if length(Tindex) > 1 && any(diff(Tindex) == 1)
        % It's a TT file
        ftype = 'tt';
        fmt1 = '%s.TT%d';
    end
    if (dashindex(end) < Tindex(end)) ...
            && (Tindex(end) < Cindex(end))
        trodenum = str2num(...
            spikename((Tindex(end)+1) : (Cindex(end)-1)) );
    end
    if isempty(trodenum)
        trodenum = 0;
    end
    OutputFileName = sprintf(fmt1, ...
        spikename( 1:(dashindex(end)-1) ), ...
        trodenum );
end

