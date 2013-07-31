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

% Last Modified by GUIDE v2.5 30-Jul-2013 13:32:28

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

%default prefs placeholders for testing:
%region.timeres = 1;
%handles.prefs = setupDetectionPreferences(region);

% Read the optional parameters
if (rem(length(varargin),2)==1)
  error('Optional parameters should always go by pairs');
else
  for i=1:2:(length(varargin)-1)
    if ~ischar (varargin{i}),
      error (['Unknown type of optional parameter name (parameter' ...
	      ' names must be strings).']);
    end
    % change the value of parameter
    switch varargin{i}
     case 'prefs'
      handles.prefs = varargin{i+1};     
     otherwise
      % Hmmm, something wrong with the parameter string
      error(['Unrecognized parameter: ''' varargin{i} '''']);
    end;
  end;
end

%handles.current = 'normal';
set(handles.popupmenu1,'String',{handles.prefs.name})

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



function prefs = getCurrentParams(hObject, eventdata, handles)
idx = get(handles.popupmenu1,'Value');  %use idx of popupmenu1 to index into {handles.prefs.name} to fetch the corresponding default params below
contents = cellstr(get(handles.popupmenu1,'String'));
prefs.name = contents{idx};

prefs.params = handles.prefs(idx).params;
for i = 1:length(prefs.params)
    switch prefs.params(i).name
        case 'hannfilterorder'
            prefs.params(i).value = str2double(get(handles.hannfilterorder , 'String'));
        case 'sd'
            prefs.params(i).value = str2double(get(handles.sd, 'String'));
        case 'sd2'
            prefs.params(i).value = str2double(get(handles.sd2, 'String'));
        case 'sd3'
            prefs.params(i).value = str2double(get(handles.sd3, 'String'));
        case 'nonfilt'
            prefs.params(i).value = str2double(get(handles.nonfilt , 'String'));
        case 'hipass'
            prefs.params(i).value = get(handles.hipass , 'String');
        case 'slidingWinStartFrame'
            prefs.params(i).value = str2double(get(handles.slidingWinStartFrame , 'String'));
        case 'block_size'
            prefs.params(i).value = str2double(get(handles.block_size , 'String'));
        case 'start_baseline'
            prefs.params(i).value = str2double(get(handles.start_baseline , 'String'));
        case 'end_baseline'
            prefs.params(i).value = str2double(get(handles.end_baseline , 'String'));
        case 'blockSizeMultiplier'
            prefs.params(i).value = str2double(get(handles.blockSizeMultiplier , 'String'));
        case 'windowAverage'
            prefs.params(i).value = get(handles.windowAverage , 'String');
        case 'baselineAverage'
            prefs.params(i).value = get(handles.baselineAverage , 'String');
        case 'maxOffsetTime'
            prefs.params(i).value = str2double(get(handles.maxOffsetTime , 'String'));
        otherwise
            error(['Unrecognized parameter: ''' prefs.params(i).name '''']);
    end
end




% --- Executes when selected object changed in unitgroup.
function unitgroup_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in unitgroup 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% %----START setup defaults ---------------------------------------------------
idx = get(hObject,'Value');  %use idx of popupmenu1 to index into {handles.prefs.name} to fetch the corresponding default params below
%--Setup default params
for i = 1:length(handles.prefs(idx).params)
    switch handles.prefs(idx).params(i).name
        case 'hannfilterorder'
            set(handles.hannfilterorder , 'String',handles.prefs(idx).params(i).value);
        case 'sd'
            set(handles.sd, 'String',handles.prefs(idx).params(i).value);
        case 'sd2'
            set(handles.sd2, 'String',handles.prefs(idx).params(i).value);
        case 'sd3'
            set(handles.sd3, 'String',handles.prefs(idx).params(i).value);
        case 'nonfilt'
            set(handles.nonfilt , 'String',handles.prefs(idx).params(i).value);
        case 'hipass'
            set(handles.hipass , 'String',handles.prefs(idx).params(i).value);
        case 'slidingWinStartFrame'
            set(handles.slidingWinStartFrame , 'String',handles.prefs(idx).params(i).value);
        case 'block_size'
            set(handles.block_size , 'String',handles.prefs(idx).params(i).value);
        case 'start_baseline'
            set(handles.start_baseline , 'String',handles.prefs(idx).params(i).value);
        case 'end_baseline'
            set(handles.end_baseline , 'String',handles.prefs(idx).params(i).value);
        case 'blockSizeMultiplier'
            set(handles.blockSizeMultiplier , 'String',handles.prefs(idx).params(i).value);
        case 'windowAverage'
            set(handles.windowAverage , 'String',handles.prefs(idx).params(i).value);
        case 'baselineAverage'
            set(handles.baselineAverage , 'String',handles.prefs(idx).params(i).value);
        case 'maxOffsetTime'
            set(handles.maxOffsetTime , 'String',handles.prefs(idx).params(i).value);
        otherwise
            error(['Unrecognized parameter: ''' handles.prefs(idx).params(i).name '''']);
    end
end
guidata(handles.figure1, handles);




% --------------------------------------------------------------------
function initialize_gui(fig_handle, handles, isreset)
% If the metricdata field is present and the reset flag is false, it means
% we are we are just re-initializing a GUI by calling it from the cmd line
% while it is up. So, bail out as we dont want to reset the data.
% if isfield(handles, 'metricdata') && ~isreset
%     return;
% end
set(handles.unitgroup, 'SelectedObject', handles.popupmenu1);
unitgroup_SelectionChangeFcn(handles.popupmenu1, [], handles)
% Update handles structure
guidata(handles.figure1, handles);



% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1
unitgroup_SelectionChangeFcn(handles.popupmenu1, [], handles)


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


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
%Export data 'handles.region' to workspace
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
prefs = getCurrentParams(hObject, eventdata, handles);
assignin('base', 'prefs', prefs)
delete(handles.figure1)



% --- Executes on button press in pushbutton3.
function pushbutton4_Callback(hObject, eventdata, handles)
%Export data 'handles.region' to workspace
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
matlabUserPath = userpath;
matlabUserPath = matlabUserPath(1:end-1);
calciumdxprefs = fullfile(matlabUserPath,'calciumdxprefs.mat');
load(calciumdxprefs)

prefs = getCurrentParams(hObject, eventdata, handles);

if strcmp(prefs.name,'custom1')  
    %--Add custom preference sets if they exist
    if exist('dxeventsPrefs','var')
        for i=1:length(dxeventsPrefs.detector)
            switch dxeventsPrefs.detector(i).name
                case 'calciumdxdettrial'
                    for j = 1:length(dxeventsPrefs.detector(i).prefs)
                        switch dxeventsPrefs.detector(i).prefs(j).name
                            case 'custom1'
                                dxeventsPrefs.detector(i).prefs(j).params = prefs.params;
                                %                             nPrefs = nPrefs+1;
                                %                             prefs(nPrefs).name = dxeventsPrefs.detector(i).prefs(j).name;
                                %                             prefs(nPrefs).params = dxeventsPrefs.detector(i).prefs(j).params;
                        end
                    end
            end
        end
    else
        dxeventsPrefs.detector(1).name = 'calciumdxdettrial';
        dxeventsPrefs.detector(1).prefs = prefs;
        
    end
end
save(calciumdxprefs, 'dxeventsPrefs','-append')











%=======uibutton group edit text default functions=========================

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
% handles.metricdata.density = density;
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
% handles.metricdata.volume = volume;
guidata(hObject,handles)




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



function hipass_Callback(hObject, eventdata, handles)
% hObject    handle to hipass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hipass as text
%        str2double(get(hObject,'String')) returns contents of hipass as a double


% --- Executes during object creation, after setting all properties.
function hipass_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hipass (see GCBO)
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
