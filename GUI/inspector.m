function varargout = inspector(varargin)
% INSPECTOR MATLAB code for inspector.fig
%      INSPECTOR, by itself, creates a new INSPECTOR or raises the existing
%      singleton*.
%
%      H = INSPECTOR returns the handle to a new INSPECTOR or the handle to
%      the existing singleton*.
%
%      INSPECTOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in INSPECTOR.M with the given input arguments.
%
%      INSPECTOR('Property','Value',...) creates a new INSPECTOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before inspector_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to inspector_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help inspector

% Last Modified by GUIDE v2.5 15-Jul-2015 10:57:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @inspector_OpeningFcn, ...
                   'gui_OutputFcn',  @inspector_OutputFcn, ...
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


% --- Executes just before inspector is made visible.
function inspector_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to inspector (see VARARGIN)

if (length (varargin) ~= 1)
    error ('GUI-Input must not contain more than the Density-Object.');
end % if

% add density-object to the gui-handles
handles.densityObj      = varargin{1};
handles.K               = size (handles.densityObj.getMu, 2);
handles.componentColors = linspecer (handles.K);
% fill up the first popupmenu with the available
set (handles.popupmenu3, 'String', ...
    cellArray2String (handles.densityObj.getRVNames()));
% disable other popupmenues
set (handles.popupmenu1, 'Enable', 'off');
set (handles.popupmenu4, 'Enable', 'off');
set (handles.slider1,    'Enable', 'off');

% Choose default command line output for inspector
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes inspector wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = inspector_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% plot density
copiedObj = handles.densityObj.cond_const (handles.RVCond, get(hObject, 'Value'));
copiedObj.marg ({handles.RV1, handles.RV2});

copiedObj.plot;


% track densities' mean and mode
hold on

mu = copiedObj.getMu;
handles.muTrack = [handles.muTrack; mu(:)'];

for k = 1:ceil (size (handles.muTrack, 2) / 2)
    ind = ((k-1) * 2 + 1);
    plot (handles.muTrack(:, ind), handles.muTrack(:, ind + 1), ...
        'Color', handles.componentColors (k, :));
end % for

mode = copiedObj.aggregation ('mode');
plot (mode(1), mode(2), 'k*');

% handles.modeTrack = [handles.modeTrack; mode(:)'];
% plot (handles.modeTrack(:, 1), handles.modeTrack(:, 2), 'k');

hold off;

% customize plot
xlim ([str2num(get(handles.edit2, 'String')), str2num(get(handles.edit3, 'String'))]);
ylim ([str2num(get(handles.edit4, 'String')), str2num(get(handles.edit5, 'String'))]);

xlabel (handles.RV1);
ylabel (handles.RV2);

set (handles.text6, 'String', sprintf ('%s = %.4f', ...
    handles.RVCond, get(hObject, 'Value')));

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1

RVNames      = cellstr (get (hObject, 'String'));
handles.RV2  = RVNames {get (hObject, 'Value')}; 

% enable the second popupmenu
set (handles.popupmenu4, 'Enable', 'on');

% fill up the first popupmenu with the still available RV-Names
RVNames (get (hObject, 'Value')) = [];
set (handles.popupmenu4, 'String', cellArray2String (RVNames));

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3

RVNames      = cellstr (get (hObject, 'String'));
handles.RV1  = RVNames {get (hObject, 'Value')}; 

% enable the second popupmenu
set (handles.popupmenu1, 'Enable', 'on');
% fill up the first popupmenu with the still available RV-Names
RVNames (get (hObject, 'Value')) = [];
set (handles.popupmenu1, 'String', cellArray2String (RVNames));

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu4.
function popupmenu4_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu4

RVNames         = cellstr (get (hObject, 'String'));
handles.RVCond  = RVNames {get (hObject, 'Value')}; 

handles.muTrack   = [];
handles.modeTrack = [];

% enable the second popupmenu
set (handles.slider1, 'Enable', 'on');

set (handles.slider1, 'Min',   str2num (get (handles.edit6, 'String')));
set (handles.slider1, 'Max',   str2num (get (handles.edit7, 'String')));
tmp = (get (handles.slider1, 'Max') - get (handles.slider1, 'Min')) / 2;
set (handles.slider1, 'Value',  tmp);

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupmenu4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- HELPER FUNCTIONS
function str = cellArray2String (varargin)
str = varargin{1};

n = numel (varargin);
for i = 2:n
    str = sprintf ('%s\n%s', str, varargin{i});
end % for



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double
xlim ([str2num(get(handles.edit2, 'String')), str2num(get(handles.edit3, 'String'))]);
ylim ([str2num(get(handles.edit4, 'String')), str2num(get(handles.edit5, 'String'))]);



% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double
xlim ([str2num(get(handles.edit2, 'String')), str2num(get(handles.edit3, 'String'))]);
ylim ([str2num(get(handles.edit4, 'String')), str2num(get(handles.edit5, 'String'))]);



% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double
xlim ([str2num(get(handles.edit2, 'String')), str2num(get(handles.edit3, 'String'))]);
ylim ([str2num(get(handles.edit4, 'String')), str2num(get(handles.edit5, 'String'))]);



% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double
xlim ([str2num(get(handles.edit2, 'String')), str2num(get(handles.edit3, 'String'))]);
ylim ([str2num(get(handles.edit4, 'String')), str2num(get(handles.edit5, 'String'))]);



% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double

set (handles.slider1, 'Min',   str2num (get (hObject, 'String')));



% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double

set (handles.slider1, 'Max',   str2num (get (hObject, 'String')));

% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
