function varargout = lfp_select(varargin)
%LFP_SELECT is a GUI that selects trials based on timestamped events in
% the trial.

%      LFP_SELECT, by itself, creates a new LFP_SELECT or raises the existing
%      singleton*.
%
%      H = LFP_SELECT returns the handle to a new LFP_SELECT or the handle to
%      the existing singleton*.
%
%      LFP_SELECT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LFP_SELECT.M with the given input arguments.
%
%      LFP_SELECT('Property','Value',...) creates a new LFP_SELECT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before lfp_select_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to lfp_select_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

%$Rev: 287 $
%$Date: 2012-10-26 20:14:45 -0400 (Fri, 26 Oct 2012) $
%$Author: dgibson $

% Edit the above text to modify the response to help lfp_select

% Last Modified by GUIDE v2.5 19-Jan-2004 11:31:02

% Begin initialization code - DO NOT EDIT
% I edited it!  I edited it!  Hih-hih, hmm-hmm-hmm, hih-hih-hih!
if ismac || isunix
    myguiname = 'lfp_selectMac';
else
    myguiname = 'lfp_select';
end
gui_Singleton = 1;
gui_State = struct('gui_Name',       myguiname, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @lfp_select_OpeningFcn, ...
                   'gui_OutputFcn',  @lfp_select_OutputFcn, ...
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



function updateTrialsDisplay(handles)
% This is an inherently perverted function because 'image' displays from
% bottom to top and then left to right as we proceed through the matrix
% elements of the image by row number and then column number (i.e. it flips
% the matrix vertically), whereas we want our trials to be displayed in
% reading order, i.e. left to right and then top to bottom.  Consequently,
% we have to construct the image in reverse order, and then transpose and
% flip it.  (There's probably a smaller canonical number of operations to
% accomplish the same effect, like 2 instead of these 3, but this works.)
lfp_declareGlobals;
numV = 100;
numH = 25;
numtot = numV * numH;
rowsPerTick = 4;
YTicks = 0 : rowsPerTick : numV;
YLabels = (numtot + 1 : -(rowsPerTick*numH) : 1);
set(handles.axes2, ...
    'GridLineStyle', '-', ...
    'YTickLabel', YLabels, ...
    'YTick', YTicks);
grid(handles.axes2, 'on');

% Calculate color code array cdata: 1 = no data, 2 = not selected, 3 =
% selected, 4 = bad trial.

% Truncate lfp_SelectedTrials if needed:
if length(lfp_SelectedTrials) > numtot
    cdata = lfp_SelectedTrials(1:numtot);
else
    cdata = lfp_SelectedTrials(1:end);
end
% Convert to color codes, mark bad trials, and reverse the array:
cdata = cdata + 2;
cdata(lfp_BadTrials) = 4;
cdata = cdata(end:-1:1);
% Pad as needed to fill image:
cdata = [ones(1, numtot - length(cdata)) ...
        cdata];
cdata = fliplr(reshape(cdata, numH, numV)');
colormap([0 0 .5; 1 0 0; 0 1 0; 0 0 0]);
% Get rid of any previous images
children = get(handles.axes2, 'Children');
for h = children'
    if strcmp(get(h, 'Type'), 'image')
        delete(h);
    end
end
image('CData', cdata, 'Parent', handles.axes2, ...
    'CDataMapping', 'direct')



% --- Executes just before lfp_select is made visible.
function lfp_select_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to lfp_select (see VARARGIN)

if ismember('initialized', fieldnames(handles))
    % do nothing
else
    % Perform initialization the first time this func is called.
    
    % Choose default command line output for lfp_select
    handles.output = hObject;
    
    % UIWAIT makes lfp_select wait for user response (see UIRESUME)
    % uiwait(handles.figure1);
    
    % Add more extra conditions
    hP = get(handles.extrafiltconj(1), 'Parent');
    offset = 5;
    for copynum = 1:2
        h = copyobj(handles.extrafiltconj(1), hP);
        handles.extrafiltconj(end + 1) = h;
        pos = get(h, 'Position');
        pos(2) = pos(2) - copynum * offset;
        set(h, 'Position', pos);
        h = copyobj(handles.extrafiltcrit(1), hP);
        handles.extrafiltcrit(end + 1) = h;
        pos = get(h, 'Position');
        pos(2) = pos(2) - copynum * offset;
        set(h, 'Position', pos);
        h = copyobj(handles.extrafiltval(1), hP);
        handles.extrafiltval(end + 1) = h;
        pos = get(h, 'Position');
        pos(2) = pos(2) - copynum * offset;
        set(h, 'Position', pos);
    end
    handles.initialized = true;
    guidata(hObject, handles);
end
updateTrialsDisplay(handles);


% --- Outputs from this function are returned to the command line.
function varargout = lfp_select_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in selectallbutton.
function selectallbutton_Callback(hObject, eventdata, handles)
% hObject    handle to selectallbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
lfp_selectByRule('true')
updateTrialsDisplay(handles);

% --- Executes on button press in unselectallbutton.
function unselectallbutton_Callback(hObject, eventdata, handles)
% hObject    handle to unselectallbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
lfp_selectByRule('false');
updateTrialsDisplay(handles);


% --- Executes on button press in trialcheckbox.
function trialcheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to trialcheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of trialcheckbox


% --- Executes during object creation, after setting all properties.
function filt1crit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filt1crit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in filt1crit.
function filt1crit_Callback(hObject, eventdata, handles)
% hObject    handle to filt1crit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns filt1crit contents as cell array
%        contents{get(hObject,'Value')} returns selected item from filt1crit


% --- Executes on button press in filterbutton.
function filterbutton_Callback(hObject, eventdata, handles)
% hObject    handle to filterbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
lfp_declareGlobals;
crit = { '' '~' };
conj = { '&&' '||' };
value = get(handles.filt1val, 'String');
condition = ['HasEvent(' value ')'];
condition = [crit{get(handles.filt1crit, 'Value')} condition];
for xnum = 1:length(handles.extrafiltconj)
    if get(handles.extrafiltconj(xnum), 'Value') == 3
        break
    else
        condition = [condition ...
                conj{get(handles.extrafiltconj(xnum), 'Value')} ];
        condition = [condition ...
                crit{get(handles.extrafiltcrit(xnum), 'Value')} ];
        value = get(handles.extrafiltval(xnum), 'String');
        condition = [condition 'HasEvent(' value ')'];
    end
end
lfp_selectByRule(condition);
updateTrialsDisplay(handles);

% --- Executes during object creation, after setting all properties.
function filt1val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filt1val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function filt1val_Callback(hObject, eventdata, handles)
% hObject    handle to filt1val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of filt1val as text
%        str2double(get(hObject,'String')) returns contents of filt1val as a double


% --- Executes during object creation, after setting all properties.
function extrafiltconj_CreateFcn(hObject, eventdata, handles)
% hObject    handle to extrafiltconj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in extrafiltconj.
function extrafiltconj_Callback(hObject, eventdata, handles)
% hObject    handle to extrafiltconj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns extrafiltconj contents as cell array
%        contents{get(hObject,'Value')} returns selected item from extrafiltconj


% --- Executes during object creation, after setting all properties.
function extrafiltcrit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to extrafiltcrit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in extrafiltcrit.
function extrafiltcrit_Callback(hObject, eventdata, handles)
% hObject    handle to extrafiltcrit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns extrafiltcrit contents as cell array
%        contents{get(hObject,'Value')} returns selected item from extrafiltcrit


% --- Executes during object creation, after setting all properties.
function extrafiltval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to extrafiltval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function extrafiltval_Callback(hObject, eventdata, handles)
% hObject    handle to extrafiltval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of extrafiltval as text
%        str2double(get(hObject,'String')) returns contents of extrafiltval as a double


