function varargout = wanderlustGUI(varargin)
% wanderlust MATLAB code for wanderlust.fig
%      wanderlust by itself, creates a new wanderlust or raises the
%      existing singleton*.
%
%      H = wanderlust returns the handle to a new wanderlust or the handle to
%      the existing singleton*.
%
%      wanderlust('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in wanderlust.M with the given input arguments.
%
%      wanderlust('Property','Value',...) creates a new wanderlust or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before wanderlust_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to wanderlust_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help wanderlust

% Last Modified by GUIDE v2.5 03-May-2012 15:59:46

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @wanderlustGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @wanderlustGUI_OutputFcn, ...
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
% --- Executes just before wanderlust is made visible.
function wanderlustGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to wanderlust (see VARARGIN)

% Choose default command line output for wanderlust
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
        set(handles.txtKNeighbors, 'String', num2str(params.k));
        set(handles.txtLNeighbors, 'String', num2str(params.l));
        set(handles.txtNumGraphs, 'String', num2str(params.num_graphs));
        set(handles.txtFlockLandmarks, 'String', num2str(params.flock_landmarks));
        set(handles.txtSNN, 'String', num2str(params.snn));
        set(handles.txtNumLandmarks, 'String', num2str(params.num_landmarks));
        metrics = get(handles.pupMetric, 'String');
        metric = find(strcmpi(metrics, params.metric));
        set(handles.pupMetric, 'Value', metric);
        voting_schemes = get(handles.pupVotingScheme, 'String');
        voting_scheme = find(strcmpi(voting_schemes, params.voting_scheme));
        set(handles.pupVotingScheme, 'Value', voting_scheme);
        normalization_schemes = get(handles.pupNormalization, 'String');
        normalization_scheme = find(strcmpi(normalization_schemes, params.normalization));
        set(handles.pupNormalization, 'Value', normalization_scheme);
        set(handles.chkBandSample, 'Value', params.band_sample);
        set(handles.chkBranch, 'Value', params.branch);
        set(handles.lstGates, 'Value', params.selected_gate);
        if (params.branch) 
            set(handles.pnlParmeter,'Title', 'Enter Number of diff-map compoenents');
            set(handles.txtNumGraphs, 'String', num2str(params.kEigs));
        end
    elseif(strcmp(varargin{i}, 'gates'))
      set(handles.lstGates, 'String', varargin{i+1});
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

% UIWAIT makes wanderlust wait for user response (see UIRESUME)
uiwait(handles.figure1);
end
% --- Outputs from this function are returned to the command line.
function varargout = wanderlustGUI_OutputFcn(hObject, eventdata, handles)
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

function chkBranch_Callback(~,~, handles)
    
    if (get(handles.chkBranch, 'Value'))
        set(handles.pnlParmeter,'Title', 'Enter Number of diff-map compoenents');
        set(handles.txtNumGraphs,'String', num2str(4));
    else
        set(handles.pnlParmeter,'Title', 'Enter Number of graphs to generate');
        set(handles.txtNumGraphs,'String', num2str(10));
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
function btnRun_Callback(hObject, eventdata, handles)

    W = runWanderlust;
        
    handles.output = W;

    % Update handles structure
    guidata(hObject, handles);

    % Use UIRESUME instead of delete because the OutputFcn needs
    % to get the updated handles structure.
    uiresume(handles.figure1);
end

function W=runWanderlust
    hgui=getappdata(0,'hwand');
    handles=guihandles(hgui);

    W.selected_gate = get(handles.lstGates, 'Value');
    W.k = str2num(get(handles.txtKNeighbors, 'String'));
    W.l = str2num(get(handles.txtLNeighbors, 'String'));
    W.num_landmarks = str2num(get(handles.txtNumLandmarks, 'String'));
    W.snn = str2num(get(handles.txtSNN, 'String'));
    W.flock_landmarks = str2num(get(handles.txtFlockLandmarks, 'String'));

    nMetric = get(handles.pupMetric, 'Value');
    strMetrics = get(handles.pupMetric, 'String');
    W.metric = strMetrics{nMetric};
        
    nVotingScheme = get(handles.pupVotingScheme, 'Value');
    strVotingSchemes = get(handles.pupVotingScheme, 'String');
    W.voting_scheme = strVotingSchemes{nVotingScheme};
        
    nNormalization = get(handles.pupNormalization, 'Value');
    strNorms = get(handles.pupNormalization, 'String');
    W.normalization = strNorms{nNormalization};
        
    W.branch = get(handles.chkBranch, 'Value');
    W.band_sample = get(handles.chkBandSample, 'Value');
    
    if (W.branch)
        W.num_graphs = 1;
    else
        W.num_graphs = str2num(get(handles.txtNumGraphs, 'String'));
    end
    W.kEigs      = str2num(get(handles.txtNumGraphs, 'String'));

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

function txtMaxDistance_Callback(hObject, eventdata, handles)
% hObject    handle to txtMaxDistance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtMaxDistance as text
%        str2double(get(hObject,'String')) returns contents of txtMaxDistance as a double
end

% --- Executes during object creation, after setting all properties.
function txtMaxDistance_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtMaxDistance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function txtKNeighbors_Callback(hObject, eventdata, handles)
% hObject    handle to txtKNeighbors (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtKNeighbors as text
%        str2double(get(hObject,'String')) returns contents of txtKNeighbors as a double
end

% --- Executes during object creation, after setting all properties.
function txtKNeighbors_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtKNeighbors (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on selection change in lstChannels.
function lstChannels_Callback(hObject, eventdata, handles)
% hObject    handle to lstChannels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns lstChannels contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lstChannels
end

% --- Executes during object creation, after setting all properties.
function lstChannels_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lstChannels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end