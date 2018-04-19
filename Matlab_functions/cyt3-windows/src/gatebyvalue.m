function varargout = gatebyvalue(varargin)
% gatebyvalue MATLAB code for gatebyvalue.fig
%      gatebyvalue by itself, creates a new gatebyvalue or raises the
%      existing singleton*.
%
%      H = gatebyvalue returns the handle to a new gatebyvalue or the handle to
%      the existing singleton*.
%
%      gatebyvalue('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in gatebyvalue.M with the given input arguments.
%
%      gatebyvalue('Property','Value',...) creates a new gatebyvalue or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gatebyvalue_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gatebyvalue_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gatebyvalue

% Last Modified by GUIDE v2.5 03-May-2012 15:59:46

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gatebyvalue_OpeningFcn, ...
                   'gui_OutputFcn',  @gatebyvalue_OutputFcn, ...
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
end
% --- Executes just before gatebyvalue is made visible.
function gatebyvalue_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gatebyvalue (see VARARGIN)

% Choose default command line output for gatebyvalue
handles.output = {};

% Update handles structure
guidata(hObject, handles);
setappdata(0,'hwand',gcf);

% Insert custom Title and Text if specified by the user
% Hint: when choosing keywords, be sure they are not easily confused 
% with existing figure properties.  See the output of set(figure) for
% a list of figure properties.

hgui=getappdata(0,'hwand');

for i=1:length(varargin)-1        
    
    if(strcmp(varargin{i}, 'title'))
      set(hObject, 'Name', varargin{i+1});
    elseif(strcmp(varargin{i}, 'params'))
        params = varargin{i+1};
        set(handles.txtbottom, 'String', num2str(params.bottom));
        set(handles.txttop, 'String', num2str(params.top));
        handles.params = params;
        % Update handles structure
        guidata(hObject, handles);
    end
end

% Determine the position of the dialog - centered on the callback figure
% if available, else, centered on the screen
FigPos=get(0,'DefaultFigurePosition');
OldUnits = get(hObject, 'Units');
set(hObject, 'Units', 'pixels');
OldPos = get(hObject,'Position');
FigWidth = OldPos(3);
FigHeight = OldPos(4);
if isempty(gcbf)
    ScreenUnits=get(0,'Units');
    set(0,'Units','pixels');
    ScreenSize=get(0,'ScreenSize');
    set(0,'Units',ScreenUnits);

    FigPos(1)=1/2*(ScreenSize(3)-FigWidth);
    FigPos(2)=2/3*(ScreenSize(4)-FigHeight);
else
    GCBFOldUnits = get(gcbf,'Units');
    set(gcbf,'Units','pixels');
    GCBFPos = get(gcbf,'Position');
    set(gcbf,'Units',GCBFOldUnits);
    FigPos(1:2) = [(GCBFPos(1) + GCBFPos(3) / 2) - FigWidth / 2, ...
                   (GCBFPos(2) + GCBFPos(4) / 2) - FigHeight / 2];
end
FigPos(3:4)=[FigWidth FigHeight];
set(hObject, 'Position', FigPos);
set(hObject, 'Units', OldUnits);

% Show a question icon from dialogicons.mat - variables questIconData
% and questIconMap
load dialogicons.mat

IconData=questIconData;
questIconMap(256,:) = get(handles.figure1, 'Color');
IconCMap=questIconMap;

set(handles.figure1, 'Colormap', IconCMap);

% Make the GUI modal
set(handles.figure1,'WindowStyle','modal')

% UIWAIT makes gatebyvalue wait for user response (see UIRESUME)
uiwait(handles.figure1);
end
% --- Outputs from this function are returned to the command line.
function varargout = gatebyvalue_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% The figure can be deleted now
delete(handles.figure1);
end

% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)

    if isequal(get(hObject, 'waitstatus'), 'waiting')
        % The GUI is still in UIWAIT, us UIRESUME
        uiresume(hObject);
    else
        handles.output = [];
        % The GUI is no longer waiting, just close it
        delete(hObject);
    end
end

% --- Executes on key press over figure1 with no controls selected.
function figure1_KeyPressFcn(hObject, eventdata, handles)

    % Check for "enter" or "escape"
    if isequal(get(hObject,'CurrentKey'),'escape')
        % User said no by hitting escape
        handles.output = [];

        % Update handles structure
        guidata(hObject, handles);

        uiresume(handles.figure1);
    end    

    if isequal(get(hObject,'CurrentKey'),'return')
        uiresume(handles.figure1);
    end   
end

% --- Executes on button press in btnRun.
function btnGate_Callback(hObject, eventdata, handles)

    hgui=getappdata(0,'hwand');
    handles=guihandles(hgui);

    output.bottom = str2num(get(handles.txtbottom, 'String'));
    output.top = str2num(get(handles.txttop, 'String'));
    handles.output = output;

    % Update handles structure
    guidata(hObject, handles);

    % Use UIRESUME instead of delete because the OutputFcn needs
    % to get the updated handles structure.
    uiresume(handles.figure1);
end

% --- Executes on button press in btnCancel.
function btnCancel_Callback(hObject, eventdata, handles)
% hObject    handle to btnCancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.output = {};

% Update handles structure
guidata(hObject, handles);

% Use UIRESUME instead of delete because the OutputFcn needs
% to get the updated handles structure.
uiresume(handles.figure1);
end


function txtbottom_Callback(hObject, eventdata, handles)
% hObject    handle to txtKNeighbors (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtKNeighbors as text
%        str2double(get(hObject,'String')) returns contents of txtKNeighbors as a double

bottom = str2double(get(hObject, 'String'));
if isnan(bottom)
    set(hObject, 'String', num2str(handles.params.bottom));
    errordlg('Input must be a number','Error');
end
end

% --- Executes during object creation, after setting all properties.
function txtbottom_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtKNeighbors (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function txttop_Callback(hObject, eventdata, handles)
% hObject    handle to txtKNeighbors (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtKNeighbors as text
%        str2double(get(hObject,'String')) returns contents of txtKNeighbors as a double

top = str2double(get(hObject, 'String'));
if isnan(top)
    set(hObject, 'String', num2str(handles.params.top));
    errordlg('Input must be a number','Error');
end
end

% --- Executes during object creation, after setting all properties.
function txttop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtKNeighbors (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

