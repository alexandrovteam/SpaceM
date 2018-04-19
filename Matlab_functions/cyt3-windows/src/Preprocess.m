function varargout = Preprocess(varargin)
% PREPROCESS MATLAB code for Preprocess.fig
%      PREPROCESS by itself, creates a new PREPROCESS or raises the
%      existing singleton*.
%
%      H = PREPROCESS returns the handle to a new PREPROCESS or the handle to
%      the existing singleton*.
%
%      PREPROCESS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PREPROCESS.M with the given input arguments.
%
%      PREPROCESS('Property','Value',...) creates a new PREPROCESS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Preprocess_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Preprocess_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

    % Edit the above text to modify the response to help Preprocess

    % Last Modified by GUIDE v2.5 21-May-2012 17:38:50

    % Begin initialization code - DO NOT EDIT
    gui_Singleton = 1;
    gui_State = struct('gui_Name',       mfilename, ...
                       'gui_Singleton',  gui_Singleton, ...
                       'gui_OpeningFcn', @Preprocess_OpeningFcn, ...
                       'gui_OutputFcn',  @Preprocess_OutputFcn, ...
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

% --- Executes just before Preprocess is made visible.
function Preprocess_OpeningFcn(hObject, eventdata, handles, varargin)
    % Choose default command line output for Preprocess
    handles.output = {1, [], 0, 0, 0, 0, ''};

    % Update handles structure
    guidata(hObject, handles);

    % Insert custom Title and Text if specified by the user
    % Hint: when choosing keywords, be sure they are not easily confused 
    % with existing figure properties.  See the output of set(figure) for
    % a list of figure properties.
    if(nargin > 3)
        for index = 1:2:(nargin-3),
            if nargin-3==index, break, end
            switch lower(varargin{index})
             case 'title'
              set(hObject, 'Name', varargin{index+1});
             case 'channelnames'
              set(handles.lstAllChannels, 'String', varargin{index+1});
              set(handles.lstAllChannels, 'Value', 3:size(varargin{index+1}, 2));
            end
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

    % Make the GUI modal
    set(handles.figure1,'WindowStyle','modal')

    % UIWAIT makes Preprocess wait for user response (see UIRESUME)
    uiwait(handles.figure1);
end

% --- Outputs from this function are returned to the command line.
function varargout = Preprocess_OutputFcn(hObject, eventdata, handles)

    % Get default command line output from handles structure
    varargout = handles.output;

    % The figure can be deleted now
    delete(handles.figure1);
end
% --- Executes on button press in btnOK.
function btnOK_Callback(hObject, eventdata, handles)

    sCofactor = get(handles.txtCofactor, 'String');
    try 
        cofactor = str2double(sCofactor);
    catch 
        handles.output = {1, [], 0, 0, 0, 0, ''};
        errordlg(springf('Cofactor must have numeric value. current value is ''%s''', sCofactor),'Invalid cofactor value');
        return;
    end
    selectedChannels = get(handles.lstAllChannels, 'Value');
    isDnagate = get(handles.chkDNAGate, 'Value');
    isSaveout= get(handles.chkSaveGate, 'Value');
    isOriginal = get(handles.chkOriginal, 'Value');
    isPrefix = get(handles.chkPrefix, 'Value');
    strPrefix = get(handles.txtPrefix, 'String');
    
    handles.output = {cofactor, selectedChannels, isDnagate,...
                      isSaveout, isOriginal, isPrefix, strPrefix};

    % Update handles structure
    guidata(hObject, handles);

    % Use UIRESUME instead of delete because the OutputFcn needs
    % to get the updated handles structure.
    uiresume(handles.figure1);
end

% --- Executes on button press in btnCancel.
function btnCancel_Callback(hObject, eventdata, handles)
    handles.output = {0, [], 0, 0, 0, 0, ''};

    % Update handles structure
    guidata(hObject, handles);

    % Use UIRESUME instead of delete because the OutputFcn needs
    % to get the updated handles structure.
    uiresume(handles.figure1);
end

% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)

    if isequal(get(hObject, 'waitstatus'), 'waiting')
        % The GUI is still in UIWAIT, us UIRESUME
        uiresume(hObject);
    else
        handles.output = {1, [], 0,0};
        % The GUI is no longer waiting, just close it
        delete(hObject);
    end
end


% --- Executes on key press over figure1 with no controls selected.
function figure1_KeyPressFcn(hObject, eventdata, handles)

    % Check for "enter" or "escape"
    if isequal(get(hObject,'CurrentKey'),'escape')
        % User said no by hitting escape
        handles.output = {1, [], 0,0};

        % Update handles structure
        guidata(hObject, handles);

        uiresume(handles.figure1);
    end    

    if isequal(get(hObject,'CurrentKey'),'return')
        uiresume(handles.figure1);
    end   
end


% --- Executes on button press in chkSaveGate.
function chkSaveGate_Callback(hObject, eventdata, handles)
gethand = handles;
isSaveGate = get(handles.chkSaveGate, 'Value');
enable = 'off';
if isSaveGate
    enable = 'on';
end
    

set(handles.chkOriginal, 'Enable', enable);
set(handles.chkPrefix, 'Enable', enable);
set(handles.txtPrefix, 'Enable', enable);
end

% --- Executes on selection change in lstAllChannels.
function lstAllChannels_Callback(hObject, eventdata, handles)
end

% --- Executes during object creation, after setting all properties.
function lstAllChannels_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end    
end

% --- Executes on selection change in lstSelectedChannels.
function lstSelectedChannels_Callback(hObject, eventdata, handles)
    ctm = get(hObject, 'UIContextMenu');
    fprintf('context menu: %d', ctm);

end

% --- Executes during object creation, after setting all properties.
function lstSelectedChannels_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
end

% --- Executes on button press in btnUnselect.
function btnUnselect_Callback(hObject, eventdata, handles)
end


function txtCofactor_Callback(hObject, eventdata, handles)
end

% --- Executes during object creation, after setting all properties.
function txtCofactor_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
end

function myCallbackFcn(hObject,jEventData,hListbox)
   % Determine the click type
   % (can similarly test for CTRL/ALT/SHIFT-click)
   if jEventData.isMetaDown  % right-click is like a Meta-button
      clickType = 'Right-click';
   else
      clickType = 'Left-click';
   end
 
   % Determine the current listbox index
   % Remember: Java index starts at 0, Matlab at 1
   mousePos = java.awt.Point(jEventData.getX, jEventData.getY);
   clickedIndex = jListbox.locationToIndex(mousePos) + 1;
   listValues = get(hListbox,'string');
   clickedValue = listValues{clickedIndex};
 
   fprintf('%s on item #%d (%s)\n', clickType, clickedIndex, clickedValue);
end  % mousePressedCallback

% --------------------------------------------------------------------
function tmp_Callback(hObject, eventdata, handles)
end

% --- Executes during object creation, after setting all properties.
function txtPrefix_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes on button press in chkOriginal.
function chkOriginal_Callback(hObject, eventdata, handles)
    set(handles.chkPrefix, 'Value', 0);
    set(handles.txtPrefix, 'Enable', 'off');    
end
% --- Executes on button press in chkPrefix.
function chkPrefix_Callback(hObject, eventdata, handles)
   set(handles.chkOriginal, 'Value', 0);
    set(handles.txtPrefix, 'Enable', 'on');
end
function txtPrefix_Callback(hObject, eventdata, handles)
end
