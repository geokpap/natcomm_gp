function varargout = lfp_getRodentTaskParams(varargin)
%taskparams = lfp_getRodentTaskParams
%   GUI for specifying rodent task parameters when saving a Rodent Clusters
%   file.  <taskparams> is a 4-element vector of values coded as in the
%   Rodent Clusters file header: [ ProcType RightCS LeftCS NoGoCS ].

% LFP_GETRODENTTASKPARAMS M-file for lfp_getRodentTaskParams.fig
%      LFP_GETRODENTTASKPARAMS, by itself, creates a new LFP_GETRODENTTASKPARAMS or raises the existing
%      singleton*.
%
%      H = LFP_GETRODENTTASKPARAMS returns the handle to a new LFP_GETRODENTTASKPARAMS or the handle to
%      the existing singleton*.
%
%      LFP_GETRODENTTASKPARAMS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LFP_GETRODENTTASKPARAMS.M with the given input arguments.
%
%      LFP_GETRODENTTASKPARAMS('Property','Value',...) creates a new LFP_GETRODENTTASKPARAMS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before lfp_getRodentTaskParams_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to lfp_getRodentTaskParams_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

%$Rev: 32 $
%$Date: 2008-12-14 16:07:41 -0500 (Sun, 14 Dec 2008) $
%$Author: dgibson $

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help lfp_getRodentTaskParams

% Last Modified by GUIDE v2.5 03-Jun-2005 15:35:46

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @lfp_getRodentTaskParams_OpeningFcn, ...
                   'gui_OutputFcn',  @lfp_getRodentTaskParams_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before lfp_getRodentTaskParams is made visible.
function lfp_getRodentTaskParams_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to lfp_getRodentTaskParams (see VARARGIN)

% Choose default command line output for lfp_getRodentTaskParams
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% Initialize
if exist('lfp_getRodentTaskParams_values.mat') == 2
    load('lfp_getRodentTaskParams_values.mat');
    set(handles.ProcTypePopup, 'Value', procType);
    set(handles.RightCSPopup, 'Value', rightCS);
    set(handles.LeftCSPopup, 'Value', leftCS);
    set(handles.NoGoCSPopup, 'Value', nogoCS);
end

% UIWAIT makes lfp_getRodentTaskParams wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = lfp_getRodentTaskParams_OutputFcn(hObject, eventdata, handles) 
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


% --- Executes on selection change in ProcTypePopup.
function ProcTypePopup_Callback(hObject, eventdata, handles)
% hObject    handle to ProcTypePopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns ProcTypePopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ProcTypePopup
if ~ismember(get(hObject, 'Value'), [2])
    % Not a No-Go task
    set(handles.NoGoCSPopup, 'Value', 4);
end

% --- Executes during object creation, after setting all properties.
function ProcTypePopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ProcTypePopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in RightCSPopup.
function RightCSPopup_Callback(hObject, eventdata, handles)
% hObject    handle to RightCSPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns RightCSPopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from RightCSPopup


% --- Executes during object creation, after setting all properties.
function RightCSPopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RightCSPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in NoGoCSPopup.
function NoGoCSPopup_Callback(hObject, eventdata, handles)
% hObject    handle to NoGoCSPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns NoGoCSPopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from NoGoCSPopup


% --- Executes during object creation, after setting all properties.
function NoGoCSPopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NoGoCSPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
error('This should not exist');


% --- Executes on selection change in LeftCSPopup.
function LeftCSPopup_Callback(hObject, eventdata, handles)
% hObject    handle to LeftCSPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns LeftCSPopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from LeftCSPopup


% --- Executes during object creation, after setting all properties.
function LeftCSPopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LeftCSPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in OKButton.
function OKButton_Callback(hObject, eventdata, handles)
% hObject    handle to OKButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
stimTypeValues = [ 1 8 0 99 ];
procType = get(handles.ProcTypePopup, 'Value');
rightCS = get(handles.RightCSPopup, 'Value');
leftCS = get(handles.LeftCSPopup, 'Value');
nogoCS = get(handles.NoGoCSPopup, 'Value');
save('lfp_getRodentTaskParams_values.mat', 'procType', 'rightCS', 'leftCS', 'nogoCS');
handles.output = [ procType
    stimTypeValues(rightCS)
    stimTypeValues(leftCS)
    stimTypeValues(nogoCS) ];
guidata(hObject, handles);
uiresume(handles.figure1);


