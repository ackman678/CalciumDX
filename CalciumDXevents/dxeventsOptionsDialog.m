function varargout = dxeventsOptionsDialog(varargin)
% DXEVENTSOPTIONSDIALOG MATLAB code for dxeventsOptionsDialog.fig
%      DXEVENTSOPTIONSDIALOG, by itself, creates a new DXEVENTSOPTIONSDIALOG or raises the existing
%      singleton*.
%
%      H = DXEVENTSOPTIONSDIALOG returns the handle to a new DXEVENTSOPTIONSDIALOG or the handle to
%      the existing singleton*.
%
%      DXEVENTSOPTIONSDIALOG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DXEVENTSOPTIONSDIALOG.M with the given input arguments.
%
%      DXEVENTSOPTIONSDIALOG('Property','Value',...) creates a new DXEVENTSOPTIONSDIALOG or raises
%      the existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before dxeventsOptionsDialog_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to dxeventsOptionsDialog_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help dxeventsOptionsDialog

% Last Modified by GUIDE v2.5 29-Jul-2013 17:41:56

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @dxeventsOptionsDialog_OpeningFcn, ...
                   'gui_OutputFcn',  @dxeventsOptionsDialog_OutputFcn, ...
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

% --- Executes just before dxeventsOptionsDialog is made visible.
function dxeventsOptionsDialog_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to dxeventsOptionsDialog (see VARARGIN)

% Choose default command line output for dxeventsOptionsDialog
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

initialize_gui(hObject, handles, false);

% UIWAIT makes dxeventsOptionsDialog wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = dxeventsOptionsDialog_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function hannfilterorder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hannfilterorder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function hannfilterorder_Callback(hObject, eventdata, handles)
% hObject    handle to hannfilterorder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hannfilterorder as text
%        str2double(get(hObject,'String')) returns contents of hannfilterorder as a double
% density = str2double(get(hObject, 'String'));
% if isnan(density)
%     set(hObject, 'String', 0);
%     errordlg('Input must be a number','Error');
% end

% Save the new hannfilterorder value
handles.metricdata.density = density;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function sd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sd_Callback(hObject, eventdata, handles)
% hObject    handle to sd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sd as text
%        str2double(get(hObject,'String')) returns contents of sd as a double
% volume = str2double(get(hObject, 'String'));
% if isnan(volume)
%     set(hObject, 'String', 0);
%     errordlg('Input must be a number','Error');
% end

% Save the new sd value
handles.metricdata.volume = volume;
guidata(hObject,handles)



% % --- Executes on button press in calculate.
% function calculate_Callback(hObject, eventdata, handles)
% % hObject    handle to calculate (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% mass = handles.metricdata.density * handles.metricdata.volume;
% set(handles.mass, 'String', mass);
% 
% % --- Executes on button press in reset.
% function reset_Callback(hObject, eventdata, handles)
% % hObject    handle to reset (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% initialize_gui(gcbf, handles, true);
% 



% --- Executes when selected object changed in unitgroup.
function unitgroup_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in unitgroup 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (hObject == handles.normal)
    set(handles.hannfilterorder, 'String', 2);
    set(handles.sd, 'String', 3);
else
    set(handles.hannfilterorder, 'String', 10);
    set(handles.sd, 'String', 3);
end

% --------------------------------------------------------------------
function initialize_gui(fig_handle, handles, isreset)
% If the metricdata field is present and the reset flag is false, it means
% we are we are just re-initializing a GUI by calling it from the cmd line
% while it is up. So, bail out as we dont want to reset the data.
if isfield(handles, 'metricdata') && ~isreset
    return;
end

% handles.metricdata.density = 0;
% handles.metricdata.volume  = 0;
% 
% set(handles.hannfilterorder, 'String', handles.metricdata.density);
% set(handles.sd,  'String', handles.metricdata.volume);
% set(handles.mass, 'String', 0);

set(handles.unitgroup, 'SelectedObject', handles.normal);

set(handles.hannfilterorder, 'String', 2);
set(handles.sd, 'String', 3);

% Update handles structure
guidata(handles.figure1, handles);



function sd2_Callback(hObject, eventdata, handles)
% hObject    handle to sd2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sd2 as text
%        str2double(get(hObject,'String')) returns contents of sd2 as a double


% --- Executes during object creation, after setting all properties.
function sd2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sd2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sd3_Callback(hObject, eventdata, handles)
% hObject    handle to sd3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sd3 as text
%        str2double(get(hObject,'String')) returns contents of sd3 as a double


% --- Executes during object creation, after setting all properties.
function sd3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sd3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nonfilt_Callback(hObject, eventdata, handles)
% hObject    handle to nonfilt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nonfilt as text
%        str2double(get(hObject,'String')) returns contents of nonfilt as a double


% --- Executes during object creation, after setting all properties.
function nonfilt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nonfilt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function highpass_Callback(hObject, eventdata, handles)
% hObject    handle to highpass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of highpass as text
%        str2double(get(hObject,'String')) returns contents of highpass as a double


% --- Executes during object creation, after setting all properties.
function highpass_CreateFcn(hObject, eventdata, handles)
% hObject    handle to highpass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function slidingWinStartFrame_Callback(hObject, eventdata, handles)
% hObject    handle to slidingWinStartFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of slidingWinStartFrame as text
%        str2double(get(hObject,'String')) returns contents of slidingWinStartFrame as a double


% --- Executes during object creation, after setting all properties.
function slidingWinStartFrame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slidingWinStartFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function blockSizeMultiplier_Callback(hObject, eventdata, handles)
% hObject    handle to blockSizeMultiplier (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of blockSizeMultiplier as text
%        str2double(get(hObject,'String')) returns contents of blockSizeMultiplier as a double


% --- Executes during object creation, after setting all properties.
function blockSizeMultiplier_CreateFcn(hObject, eventdata, handles)
% hObject    handle to blockSizeMultiplier (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function block_size_Callback(hObject, eventdata, handles)
% hObject    handle to block_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of block_size as text
%        str2double(get(hObject,'String')) returns contents of block_size as a double


% --- Executes during object creation, after setting all properties.
function block_size_CreateFcn(hObject, eventdata, handles)
% hObject    handle to block_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function start_baseline_Callback(hObject, eventdata, handles)
% hObject    handle to start_baseline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of start_baseline as text
%        str2double(get(hObject,'String')) returns contents of start_baseline as a double


% --- Executes during object creation, after setting all properties.
function start_baseline_CreateFcn(hObject, eventdata, handles)
% hObject    handle to start_baseline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function end_baseline_Callback(hObject, eventdata, handles)
% hObject    handle to end_baseline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of end_baseline as text
%        str2double(get(hObject,'String')) returns contents of end_baseline as a double


% --- Executes during object creation, after setting all properties.
function end_baseline_CreateFcn(hObject, eventdata, handles)
% hObject    handle to end_baseline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maxOffsetTime_Callback(hObject, eventdata, handles)
% hObject    handle to maxOffsetTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxOffsetTime as text
%        str2double(get(hObject,'String')) returns contents of maxOffsetTime as a double


% --- Executes during object creation, after setting all properties.
function maxOffsetTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxOffsetTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function windowAverage_Callback(hObject, eventdata, handles)
% hObject    handle to windowAverage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of windowAverage as text
%        str2double(get(hObject,'String')) returns contents of windowAverage as a double


% --- Executes during object creation, after setting all properties.
function windowAverage_CreateFcn(hObject, eventdata, handles)
% hObject    handle to windowAverage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function baselineAverage_Callback(hObject, eventdata, handles)
% hObject    handle to baselineAverage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of baselineAverage as text
%        str2double(get(hObject,'String')) returns contents of baselineAverage as a double


% --- Executes during object creation, after setting all properties.
function baselineAverage_CreateFcn(hObject, eventdata, handles)
% hObject    handle to baselineAverage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
