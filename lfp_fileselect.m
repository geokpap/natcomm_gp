function varargout = lfp_fileselect(varargin)
%LFP_FILESELECT M-file for lfp_fileselect.fig
% selectedfiles = lfp_fileselect
% selectedfiles = lfp_fileselect(filetypes, message)
% <filetypes> may be:
%   'all': default
%   'csc': CSC files only
%   'evt': event files only
% <message> is displayed when soliciting a specific filetype(s), i.e. when
% <filetypes> is not 'all'.  In this case the background of the GUI also
% assumes a filetype-specific color.
% Returns cell column vector of relative filenames including extension.
% Sets lfp_DataDir as a side-effect.
%
% To determine selected file type after returning from GUI, examine
% lfp_fileselect_FileTypesPopup_string in lfp_fileselect.mat.
% lfp_fileselect_FileTypesPopup_value is also stored there, but the
% meanings of the values may change as new types are added, so its value
% should normally only be used in this M-file.

%      LFP_FILESELECT, by itself, creates a new LFP_FILESELECT or raises the existing
%      singleton*.
%
%      H = LFP_FILESELECT returns the handle to a new LFP_FILESELECT or the handle to
%      the existing singleton*.
%
%      LFP_FILESELECT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LFP_FILESELECT.M with the given input arguments.
%
%      LFP_FILESELECT('Property','Value',...) creates a new LFP_FILESELECT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before lfp_fileselect_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to lfp_fileselect_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

%$Rev: 174 $
%$Date: 2010-09-28 20:53:41 -0400 (Tue, 28 Sep 2010) $
%$Author: dgibson $

% Edit the above text to modify the response to help lfp_fileselect

% Last Modified by GUIDE v2.5 11-Mar-2004 11:59:23

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @lfp_fileselect_OpeningFcn, ...
                   'gui_OutputFcn',  @lfp_fileselect_OutputFcn, ...
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


% --- Executes just before lfp_fileselect is made visible.
function lfp_fileselect_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to lfp_fileselect (see VARARGIN)

% Choose default command line output for lfp_fileselect
handles.output = [];

if length(varargin) == 0
    set(handles.MessageLabel, 'String', '');
    handles.showfiletypes = 'all';
elseif length(varargin) == 2
    set(handles.MessageLabel, 'String', varargin{2});
    handles.showfiletypes = varargin{1};
else
    error('lfp_fileselect:badarglist', ...
        'There must be zero or two arguments');
end

% Update handles structure
guidata(hObject, handles);

switch handles.showfiletypes
    % IMPORTANT: values here depend on lfp_fileselect_FileTypesPopup_value,
    % must be modified when re-ordering
    % lfp_fileselect_FileTypesPopup_string
    case 'all'
    case 'csc'
        filetypestrings = get(handles.FileTypesPopup, 'String');
        filetypestrings = filetypestrings([1 2 3]);
        set(handles.FileTypesPopup, 'String', filetypestrings);
        set(handles.figure1, 'Color', [.6 .6 .8]);
    case 'evt'
        filetypestrings = get(handles.FileTypesPopup, 'String');
        filetypestrings = filetypestrings([4 5 6 7]);
        set(handles.FileTypesPopup, 'String', filetypestrings);
        set(handles.FileList, 'Max', 1);
        set(handles.figure1, 'Color', [.6 .8 .6]);
end
loadGUIState(handles);
handles = selectDefaultFiletype(handles);
updateFileList(handles);

% UIWAIT makes lfp_fileselect wait for user response (see UIRESUME)
% The really annoying thing about this is that it has to be commented-out
% in order to debug, but then it has to execute in order for the figure not
% to automatically close itself.
uiwait(handles.figure1);


function handles = selectDefaultFiletype(handles)
% IMPORTANT:  This function changes the state of <handles> and saves the
% changed state in the GUI, thus invalidating the old value of <handles>
% that was passed in.  As a convenience, it also returns the new value so
% that it is not necessary to call guidata again.  Note also that this must
% be called AFTER loadGUIState in order to have any effect.
switch handles.showfiletypes
    % IMPORTANT: values here depend on the "switch handles.showfiletypes"
    % in lfp_fileselect_OpeningFcn
    case 'all'
    case 'csc'
    case 'evt'
        % If there are .EVTSAV files, then this should be the default type
        lfp_declareGlobals;
        files = dir(lfp_DataDir);
        exts = {};
        for fidx = 1:length(files)
            [pathstr, name, exts{fidx}] = fileparts(files(fidx).name);
        end
        if any(ismember(upper(exts), '.EVTSAV'))
            set(handles.FileTypesPopup, 'Value', 4);
        elseif any(ismember(upper(exts), '.NEV'))
            set(handles.FileTypesPopup, 'Value', 2);
        elseif any(ismember(upper(exts), '.DAT'))
            set(handles.FileTypesPopup, 'Value', 3);
        else
            set(handles.FileTypesPopup, 'Value', 1);
        end
        % Update handles structure
        guidata(handles.FileTypesPopup, handles);
end


function updateFileList(handles)
% Applies the current file filter to all files in lfp_DataDir
lfp_declareGlobals;
filenamesIncSize = 256;
filenames = cell(filenamesIncSize, 1);
numfiles = 0;
set(handles.DirLabel, 'String', lfp_DataDir);
files = dir(lfp_DataDir);
value = get(handles.FileTypesPopup, 'Value');
strings = get(handles.FileTypesPopup, 'String');
for fnum = (1:length(files))
    if (~ files(fnum).isdir)
        [pathstr,name,ext] = fileparts(files(fnum).name);
        listThisFile = false;
        switch strings{value}
            case 'Events (saved - *.EVTSAV)'
                if strcmpi(ext, '.EVTSAV') listThisFile = true;
                end
            case 'CSC (new - *.NCS)'
                if strcmpi(ext, '.NCS') listThisFile = true;
                end
            case {'CSC (old - *.DAT)' 'Events (old - *.DAT)'}
                if strcmpi(ext, '.DAT') listThisFile = true;
                end
            case {'CSC (Matlab - *.MAT)' 'Events (Matlab - *.MAT)'}
                if strcmpi(ext, '.MAT') listThisFile = true;
                end
            case 'Single Cut Cluster (*.T)'
                if strcmpi(ext, '.T') listThisFile = true;
                end
            case 'Single Cut Cluster (*.MAT)'
                if strcmpi(ext, '.MAT') listThisFile = true;
                end
            case 'Multiple Cut Cluster (*.MAT)'
                if strcmpi(ext, '.MAT') listThisFile = true;
                end
            case 'Naotaka Clusters (*.DWD)'
                if strcmpi(ext, '.DWD') listThisFile = true;
                end
            case 'Neuralynx Video (*.NVT, *.DAT)'
                if strcmpi(ext, '.NVT') || strcmpi(ext, '.DAT')
                    listThisFile = true;
                end
            case 'Rodent Clusters (*.Tnn, *.TTn)'
                if ~isempty(regexpi(ext, '^\.TT[0-9]$')) ...
                        || ~isempty(regexpi(ext, '^\.T[0-9][0-9]$'))
                    listThisFile = true;
                end
            case 'Rodent Tracker (*.DAT)'
                if strcmpi(ext, '.DAT')
                    listThisFile = true;
                end
            case 'Single Electrode (*.NSE)'
                if strcmpi(ext, '.NSE') listThisFile = true;
                end
            case 'Events (new - *.NEV)'
                if strcmpi(ext, '.NEV') listThisFile = true;
                end
        end
        if listThisFile
            numfiles = numfiles + 1;
            if numfiles > length(filenames)
                filenames = [ filenames cell(filenamesIncSize, 1) ];
            end
            filenames(numfiles) = { files(fnum).name }; 
        end
    end
end
set(handles.FileList, 'Value', 1);
set(handles.FileList, 'String', filenames(1:numfiles));


function saveGUIState(handles)
if exist('lfp_fileselect.mat') == 2
    pathname = which('lfp_fileselect.mat');
else
    pathname = 'lfp_fileselect.mat';
end
lfp_fileselect_FileTypesPopup_value = get(handles.FileTypesPopup, 'Value');
strings = get(handles.FileTypesPopup, 'String');
lfp_fileselect_FileTypesPopup_string = ...
    strings(lfp_fileselect_FileTypesPopup_value);
save(pathname, 'lfp_fileselect_FileTypesPopup_value', ...
    'lfp_fileselect_FileTypesPopup_string');


function loadGUIState(handles)
if exist('lfp_fileselect.mat') == 2
    load('lfp_fileselect.mat');
    if lfp_fileselect_FileTypesPopup_value <= length(get(handles.FileTypesPopup, 'String'))
        set(handles.FileTypesPopup, ...
            'Value', lfp_fileselect_FileTypesPopup_value );
    else
        set(handles.FileTypesPopup, 'Value', 1);
    end
end


% --- Outputs from this function are returned to the command line.
function varargout = lfp_fileselect_OutputFcn(hObject, eventdata, handles)
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
function FileList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FileList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in FileList.
function FileList_Callback(hObject, eventdata, handles)
% hObject    handle to FileList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns FileList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from FileList


% --- Executes on button press in OpenFileButton.
function OpenFileButton_Callback(hObject, eventdata, handles)
% hObject    handle to OpenFileButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
files = get(handles.FileList, 'String');
value = get(handles.FileList, 'Value');
if value <= length(files)
    handles.output = files(value);
    guidata(hObject, handles);
else
    CancelButton_Callback(hObject, eventdata, handles);
    return
end
saveGUIState(handles);
uiresume(handles.figure1);


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
    handles = selectDefaultFiletype(handles);
    updateFileList(handles);
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
updateFileList(handles);
saveGUIState(handles);

