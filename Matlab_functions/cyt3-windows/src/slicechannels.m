function varargout = slicechannels(varargin)
% SLICECHANNELS MATLAB code for slicechannels.fig
%      SLICECHANNELS, by itself, creates a new SLICECHANNELS or raises the existing
%      singleton*.
%
%      H = SLICECHANNELS returns the handle to a new SLICECHANNELS or the handle to
%      the existing singleton*.
%
%      SLICECHANNELS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SLICECHANNELS.M with the given input arguments.
%
%      SLICECHANNELS('Property','Value',...) creates a new SLICECHANNELS or raises
%      the existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before slicechannels_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to slicechannels_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help slicechannels

% Last Modified by GUIDE v2.5 22-Apr-2014 13:02:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @slicechannels_OpeningFcn, ...
                   'gui_OutputFcn',  @slicechannels_OutputFcn, ...
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

% --- Executes just before slicechannels is made visible.
function slicechannels_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to slicechannels (see VARARGIN)

% Choose default command line output for slicechannels
handles.output = [];

% Update handles structure
guidata(hObject, handles);

initialize_gui(hObject, handles, false);
handles.output = [];

% UIWAIT makes slicechannels wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = slicechannels_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
% The figure can be deleted now
delete(handles.figure1);


% --- Executes during object creation, after setting all properties.
function txtNumSlices_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtNumSlices (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtNumSlices_Callback(hObject, eventdata, handles)
% hObject    handle to txtNumSlices (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtNumSlices as text
%        str2double(get(hObject,'String')) returns contents of txtNumSlices as a double
num_slices = str2num(get(hObject, 'String'));
if isnan(num_slices)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% Save the new txtNumSlices value
handles.slice.num_slices = num_slices;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function txtOverlap_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtOverlap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtOverlap_Callback(hObject, eventdata, handles)
% hObject    handle to txtOverlap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtOverlap as text
%        str2double(get(hObject,'String')) returns contents of txtOverlap as a double
overlap = str2double(get(hObject, 'String'));
if isnan(overlap) || overlap < .5 || overlap > handles.slice.num_slices-1
    set(hObject, 'String', 0);
    errordlg('Input must be a number between .5 to number of (slices-1)','Error');
end

% Save the new txtOverlap value
handles.slice.overlap = overlap;
guidata(hObject,handles)

% --- Executes on button press in calculate.
function calculate_Callback(hObject, eventdata, handles)
% hObject    handle to calculate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    
    handles.output = handles.slice;

    % Update handles structure
    guidata(hObject, handles);

    % Use UIRESUME instead of delete because the OutputFcn needs
    % to get the updated handles structure.
    uiresume(handles.figure1);

% --- Executes on button press in reset.
function reset_Callback(hObject, eventdata, handles)
% hObject    handle to reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    
    handles.output = [];

    % Update handles structure
    guidata(hObject, handles);

    % Use UIRESUME instead of delete because the OutputFcn needs
    % to get the updated handles structure.
    uiresume(handles.figure1);

% --- Executes when selected object changed in unitgroup.
function unitgroup_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in unitgroup 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function initialize_gui(fig_handle, handles, isreset)
% If the metricdata field is present and the reset flag is false, it means
% we are we are just re-initializing a GUI by calling it from the cmd line
% while it is up. So, bail out as we dont want to reset the data.
if isfield(handles, 'metricdata') && ~isreset
    return;
end

handles.slice.num_slices = 10;
handles.slice.overlap  = 2;

set(handles.txtNumSlices, 'String', handles.slice.num_slices);
set(handles.txtOverlap,  'String', handles.slice.overlap);

handles.output = handles.slice;
% Update handles structure
guidata(handles.figure1, handles);
