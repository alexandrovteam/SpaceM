function akZoom(handles, mButtons_ZoomPanReset)
% allows direct zooming and panning with the mouse in 2D plots.
%
%  Default mouse button functions (can be customized, see below):
%         Scroll wheel:  zoom in/out
%    Left Mouse Button:  select an ROI rectangle to zoom in
%  Middle Mouse Button:  pan view
%   Right Mouse Button:  reset view to default view
%
% SYNTAX:
%   akZoom
%   akZoom(handles)
%   akZoom(handles, mButtons_ZoomPanReset)
%   akZoom(mButtons_ZoomPanReset)
%
% DESCRIPTION:
%   akZoom() activates mouse control for all axes-objects in the current figure.
%
%   akZoom(handles) activates mouse control for the axes given by handles.
%     handles can be:
%     a) a single axis handle.
%     b) an array of axis handles. In this case all axes will be linked,
%        i.e. panning and zooming in one of the axes will affect the others as well.
%     c) a cell array consisting of axis handles and/or arrays of axis
%        handles. In this case mouse control is activated for all axes inside
%        the cell but linking is only activated for axes in a common array.
%     d) one of the following strings:
%        i)  'all': activates mouse control for all axes in all open figures.
%            The axes will not be linked except for axes belonging to the same
%            plotyy plot.
%        ii) 'all_linked': activates mouse control for all axes in all open
%            figures and links all axes.
%
%   akZoom(handles, mButtons_ZoomPanReset) and akZoom(mButtons_ZoomPanReset)
%     will use akZoom with customized mouse buttons. 
%     mButtons_ZoomPanReset can be:
%     'lmr', 'lrm', 'mlr', 'mrl', 'rlm' or 'rml'
%     where the letters stand for the buttons left, middle and right. 
%     So 'lmr' means left=ROI-zoom, middle=pan and right=reset.
%     Note: If you want to use a certain mouse button pattern as default,
%     just change mButtons_ZoomPanReset_default in the "Settings" section
%     below.
%
%
% EXAMPLES:
% %%  Simple Plot
%     figure
%     plot(10:24,rand(1,15));
%     akZoom();
% 
% %%  Subplots (independent axes)
%     figure
%     for k = 1:4
%       y = rand(1,15);
%       subplot(2, 2, k);
%       plot(y);
%     end
%     akZoom();
%
% %%  Subplots (mixture of linked and indipendent axes)
%     figure
%     ax = NaN(4,1);
%     for k = 1:4
%       y = rand(1,15);
%       ax(k) = subplot(2, 2, k);
%       plot(y);
%     end
%     akZoom({[ax(1),ax(3)],ax(2),[ax(3),ax(4)]});
%
% %%  Different figures (linked)
%     figure;
%     imshow(imread('peppers.png'));
%     ax1 = gca;
%     figure
%     imshow(rgb2gray(imread('peppers.png')));
%     ax2 = gca;
%     akZoom([ax1,ax2]);
%
% %% Find more examples in the file "akZoom_examples.m"
%
%
% KNOWN BUGS:
%   Axes linked manually with linkaxes() are unlinked by akZoom
%   Solve this in the future by: 
%   h = getappdata(h_ax, 'graphics_linkaxes');
%   linked_axes = h.Targets; ...
%
%
% Author: Alexander Kessel
% Affiliation: Max-Planck-Institut für Quantenoptik, Garching, München
% Contact : alexander.kessel <at> mpq.mpg.de
% Revision: 02. August 2015
%
%
% CHANGELOG:
%   2016-09-14 fix bug for fast scrolling with mouse wheel
%
%
% CREDITS: 
% - Rody P.S. Oldenhuis for his mouse_figure function which served as the
%   template for akZoom
% - Benoit Botton for reporting on an issue with axes nested in GUIs
% - Anne-Sophie Girard-Guichon for reporting and fixing a bug occuring when 
%   fast scrolling with mouse wheel

%% Settings
turnOff_PlotyyPositionListener = true; % default: true
% Turn off the PlotyyPositionListener which Matlab creates together with the
% plotyy axes. Usually it updates the position of a plotyy-axis if the other
% axis has been changed. However, this callback sometimes (for older Matlab
% versions?) raises errors in combination with the mouse-pan/zoom of akZoom.
% The drawback of turning off the PlotyyPositionListener is that when
% resizing a GUI window, the plotyy axes in this window do not resize properly
% anymore.
% Since the use of plotyy with older Matlab versions is more likely than
% the use of plotyy in GUI windows, the default value is true. 

mButtons_ZoomPanReset_default = 'lmr'; % default: 'lmr'
% Default mapping of mouse buttons for zoom, pan and reset (in this order)
% possible values: 'lmr','lrm','mlr','mrl','rlm','rml'

wheel_zoomFactor = 20; % default: 20
% zoom factor per wheel tick in [%] -> higher value means faster zooming

%% Parse input arguments
if nargin == 0
  % no input argument
  [h_ax, h_ax_linking, h_fig] = parse_handles([]);
  [zoom_button, pan_button, reset_button] = parse_mButtons_ZoomPanReset(mButtons_ZoomPanReset_default);
  
elseif nargin == 1
  if ischar(handles) && sum(ismember({'lmr','lrm','mlr','mrl','rlm','rml'}, handles))
    % input argument = mButtons_ZoomPanReset
    mButtons_ZoomPanReset = handles;
    [h_ax, h_ax_linking, h_fig] = parse_handles([]);
    [zoom_button, pan_button, reset_button] = parse_mButtons_ZoomPanReset(mButtons_ZoomPanReset);
  else
    % input argument = axis handles
    [h_ax, h_ax_linking, h_fig] = parse_handles(handles);
    [zoom_button, pan_button, reset_button] = parse_mButtons_ZoomPanReset(mButtons_ZoomPanReset_default);
  end
  
elseif nargin == 2
  % input arguments = axis handles, mButtons_ZoomPanReset
  [h_ax, h_ax_linking, h_fig] = parse_handles(handles);
  [zoom_button, pan_button, reset_button] = parse_mButtons_ZoomPanReset(mButtons_ZoomPanReset);
  
else
  error('akZoom:invalid_number_of_arguments', 'Invalid number of input arguments. Check description of function akZoom.');
end

%% Perform initialization

% Turn off matlab-linking for all axes. akZoom is linking axes on its own.
linkaxes(h_ax, 'off');
if turnOff_PlotyyPositionListener % see section "Settings" for details
  for j = 1:numel(h_ax)
    if isfield(getappdata(h_ax(j)), 'graphicsPlotyyPositionListener');
      rmappdata(h_ax(j),'graphicsPlotyyPositionListener');
    end
  end
end

% Initialize variables for use across all nested functions
cx = []; % clicked x-coordinate
cy = []; % clicked y-coordinate
mode = ''; % mouse control mode, e.g. 'pan'
ROI = []; % This will later hold the rectangle object that marks the zoom area.

% save original limits
original_xlim = NaN(numel(h_ax),2);
original_ylim = NaN(size(original_xlim));
for j=1:numel(h_ax)
  original_xlim(j,:) = get(h_ax(j), 'xlim');
  original_ylim(j,:) = get(h_ax(j), 'ylim');
end

tX = timer('StartDelay',0.5);
tY = timer('StartDelay',0.5);

% set callbacks for all figures
for j=1:numel(h_fig)
  % switch off zoom and pan modes in case they were on
  zoom(h_fig(j), 'off')
  pan(h_fig(j), 'off')
  % set callbacks
  set(h_fig(j), ...
    'WindowScrollWheelFcn' , @scroll_zoom,...
    'WindowButtonDownFcn'  , @MouseDown,...
    'WindowButtonUpFcn'    , @MouseUp,...
    'WindowButtonMotionFcn', @MouseMotion);
end

%% Mouse callback functions
  function scroll_zoom(varargin)
    h = hittest(gcf);
    if isempty(h), return, end
    try
      switch get(h,'Type')
        case 'axes'
          currAx = h;
        case 'image'
          currAx = get(h,'Parent');
        case 'line'
          currAx = get(h,'Parent');
        case 'patch'
          currAx = get(h,'Parent');
        otherwise
          return
      end
    catch
      % in Matlab2015 and newer an opaque object can be returned by the
      % hittest, which results in an error when trying to apply the 'get'
      % method on it
      return
    end
    if ~any(currAx == h_ax), return, end
    [x,y] = getAbsCoords(currAx);
    if ~coordsWithinLimits(currAx,x,y), return, end
    [x_rel, y_rel] = abs2relCoords(currAx, x, y);
    sc = varargin{2}.VerticalScrollCount;
    zoomFactor = (abs(sc)*(1+wheel_zoomFactor/100))^sign(sc);
    if zoomFactor ~= 0 % could happen when fast scrolling 
      for i = affectedAxes(currAx)
        new_xlim_rel = ([0,1] - x_rel) * zoomFactor + x_rel;
        new_ylim_rel = ([0,1] - y_rel) * zoomFactor + y_rel;
        [new_xlim(1), new_ylim(1)] = rel2absCoords(h_ax(i), new_xlim_rel(1), new_ylim_rel(1));
        [new_xlim(2), new_ylim(2)] = rel2absCoords(h_ax(i), new_xlim_rel(2), new_ylim_rel(2));
        setNewLimits(h_ax(i), new_xlim, new_ylim)
      end
    end
  end

  function MouseDown(varargin)
    currAx = h_ax( get(gcf, 'currentaxes') == h_ax);
    if isempty(currAx), return, end %return if the current axis is not one of the axes in h_ax
    if cursorOverOtherObject, return, end %return if the cursor is above any other object (e.g. a legend)
    if strcmp(mode, 'selectROI') %if we are in zoom-mode and another butten is clicked, delete ROI
      delete(ROI);
    end
    mode = '';
    [x, y] = getAbsCoords(currAx);
    if ~coordsWithinLimits(currAx,x,y), return, end
    % save clicked coordinates to cx and cy
    cx = x;   cy = y;
    switch lower(get(gcf, 'selectiontype'))
      case pan_button;
        mode = 'pan';
      case reset_button
        for i = affectedAxes(currAx)
          set(h_ax(i), 'Xlim', original_xlim(i,:), 'Ylim', original_ylim(i,:));
        end
      case zoom_button
        mode = 'selectROI';
        % create ROI rectangle object
        w = realmin; % must be >0
        h = realmin; % must be >0
        ROI = rectangle('Position',[cx,cy,w,h]);
    end
  end

  function MouseUp(varargin)
    currAx = h_ax( get(gcf, 'currentaxes') == h_ax);
    if isempty(currAx), return, end %return if the current axis is not one of the axes specified in h_ax
    if strcmp(mode, 'selectROI')
      % get corner points of ROI in relative coordinates
      pos = get(ROI,'Position');
      if pos(3)>realmin && pos(4)>realmin %check if ROI has valid size
        x = [pos(1), pos(1)+pos(3)];
        y = [pos(2), pos(2)+pos(4)];
        [x_rel1, y_rel1] = abs2relCoords(currAx, x(1), y(1));
        [x_rel2, y_rel2] = abs2relCoords(currAx, x(2), y(2));
        for i = affectedAxes(currAx)
          % calc absolute coordinates of ROI corners
          [x1, y1] = rel2absCoords(h_ax(i), x_rel1, y_rel1);
          [x2, y2] = rel2absCoords(h_ax(i), x_rel2, y_rel2);
          new_xlim = sort([x1, x2]);
          new_ylim = sort([y1, y2]);
          setNewLimits(h_ax(i), new_xlim, new_ylim)
        end
      end
      delete(ROI);
    end
    mode = '';
    cx = [];
    cy = [];
    set(gcf,'selected','on')
  end

  function MouseMotion(varargin)
    if isempty(cx), return, end % return if there is no clicked point set
    currAx = h_ax( get(gcf, 'currentaxes') == h_ax);
    if isempty(currAx), return, end %return if the current axis is not one of the axes specified in h_ax
    [x,y] = getAbsCoords(currAx);
    if strcmp(mode, 'pan')
      if ~coordsWithinLimits(currAx,x,y), return, end %return if we are outside of limits
      % calc change in position
      [x_rel, y_rel] = abs2relCoords(currAx, x, y);
      [cx_rel, cy_rel] = abs2relCoords(currAx, cx, cy);
      delta_x_rel = x_rel - cx_rel;
      delta_y_rel = y_rel - cy_rel;
      for i = affectedAxes(currAx)
        % set new limits
        [new_xlim(1), new_ylim(1)] = rel2absCoords(h_ax(i), -delta_x_rel, -delta_y_rel);
        [new_xlim(2), new_ylim(2)] = rel2absCoords(h_ax(i), 1-delta_x_rel, 1-delta_y_rel);
        setNewLimits(h_ax(i), new_xlim, new_ylim);
      end
      % save new position
      [cx,cy] = getAbsCoords(currAx);
    elseif strcmp(mode, 'selectROI')
      % if mouse is ouside limits, adjust coords
      xlim = get(currAx, 'xlim');
      ylim = get(currAx, 'ylim');
      x = max(x,xlim(1));
      x = min(x,xlim(2));
      y = max(y,ylim(1));
      y = min(y,ylim(2));
      % resize ROI rectangle
      w = max(abs(x-cx), realmin); % must be >0
      h = max(abs(y-cy), realmin); % must be >0
      set(ROI, 'Position', [min(cx,x), min(cy,y), w, h]);
    else % no mode
      return;
    end
  end



%% supporting functions
  function [zoom_button, pan_button, reset_button] = parse_mButtons_ZoomPanReset(mButtons_ZoomPanReset)
    % In this function we map the mouse buttons
    if sum(ismember({'lmr','lrm','mlr','mrl','rlm','rml'}, mButtons_ZoomPanReset))
      l = 'normal';
      m = 'extend';
      r = 'alt';
      eval(['zoom_button = ' mButtons_ZoomPanReset(1) ';']);
      eval(['pan_button = ' mButtons_ZoomPanReset(2) ';']);
      eval(['reset_button = ' mButtons_ZoomPanReset(3) ';']);
    else
      error('akZoom:invalid_mButtons_ZoomPanReset_input', ['Invalid mButtons_ZoomPanReset input: "' mButtons_ZoomPanReset '"']);
    end
  end

  function [h_ax, h_ax_linking, h_fig] = parse_handles(handles)
    % In this function we provide the three crucial variables:
    % h_ax is an array of unique axes handles: [h1, h2, h3, ...]
    % h_ax_linking is a cell array containing axes handles and/or axes handle arrays
    %   as elements. Axes handles within one axes handle array are later treated
    %   as linked, i.e. panning/zooming one axis affects the other axes in the
    %   array as well: {h1, [h2,h3], [h2,h4], [h5,h6,h7], h8, ...}
    % h_fig is an array of unique figure handles specifying the parents of the axes in h_ax
    
    % generate h_ax
    if isempty(handles)
      % find all axes in current figure
      h_fig = get(0,'CurrentFigure');
      if isempty(h_fig)
        error('akZoom:no_open_figure', 'There is no open figure.');
      end
      h_ax = findall(h_fig,'type','axes','-not','Tag','legend','-not','Tag','Colorbar');
      h_ax_linking = 'all_unlinked';
    elseif ischar(handles)
      if strcmp(handles,'all')
        h_ax = findall(0,'type','axes','-not','Tag','legend','-not','Tag','Colorbar');
        h_ax_linking = 'all_unlinked';
      elseif strcmp(handles,'all_linked')
        h_ax = findall(0,'type','axes','-not','Tag','legend','-not','Tag','Colorbar');
        h_ax_linking = 'all_linked';
      else
        error('akZoom:invalid_input', ['Invalid input string: "' handles '". Check description of function akZoom.']);
      end
    elseif iscell(handles) % cell array with linking information
      h_ax_linking = handles;
      h_ax = unique(cat(2,handles{:}));
    else % handle array -> all axes linked
      h_ax = unique(handles);
      h_ax_linking = 'all_linked';
    end
    % check h_ax
    if isempty(h_ax)
      error('akZoom:no_axis_found','No axis found.')
    elseif any(~ishandle(h_ax))
      error('akZoom:is_not_a_handle','One of the axis handles is not a valid handle.');
    elseif any(~is2D(h_ax))
      error('akZoom:axis_not_2D', 'At least one of the axes is not a 2D-axis. akZoom does not support 3D-axes so far.');
    end
    
    % generate h_ax_linking
    if ischar(h_ax_linking)
      if strcmp(h_ax_linking, 'all_linked')
        h_ax_linking = {h_ax}; % all axes linked
      elseif strcmp(h_ax_linking, 'all_unlinked')
        h_ax_linking = linkPlotyyAxes(h_ax); % all axes unlinked except plotyy-axes
      else
        error('akZoom:invalid_string', ['Invalid string for h_ax_linking: ' h_ax_linking]);
      end
    else
      % do nothing: h_ax_linking was already specified
    end
    
    % generate h_fig
    h_fig = NaN(1, numel(h_ax));
    for i=1:numel(h_ax)
      % figures are usually the parent objects of the axes
      fig = get(h_ax(i), 'Parent');
      % in nested GUIs the figure object might be on a higher inheritance level 
      while ~strcmp(get(fig,'type'), 'figure')
        fig = get(fig,'parent');
        if isempty(fig)
          error('akZoom:no_parent_figure',['No parent figure found for axis ' num2str(h_ax(i))]);
        end
      end
      h_fig(i) = fig;
    end
    if iscell(h_fig)
      h_fig = unique([h_fig{:}]);
    end
  end

  function [x, y, z] = getAbsCoords(h_ax)
    crd = get(h_ax, 'CurrentPoint');
    x = crd(2,1);
    y = crd(2,2);
    z = crd(2,3);
  end

  function tf = cursorOverOtherObject()
    tf = false;
    % check if cursor is over legend
    ax = overobj('axes');
    if ~isempty(ax) && strcmp(get(ax,'Tag'),'legend')
        tf = true;
    end
    % check if cursor is over DataTip object
    h_dt = findall(h_fig,'Tag','DataTipMarker');
    if ~isempty(h_dt)
      for i=1:numel(h_fig)
        h = hittest(h_fig(i));
        if any(h == h_dt)
          tf = true;
          return
        end
      end
    end
  end

  function [x_rel, y_rel] = abs2relCoords(h_ax, x, y)
    XLim = get(h_ax, 'xlim');
    if strcmp(get(h_ax, 'XScale'), 'log')
      x_rel = ( log(x) - log(XLim(1)) ) / ( log(XLim(2)) - log(XLim(1)) );
    else
      x_rel = (x-XLim(1))/(XLim(2)-XLim(1));
    end
    YLim = get(h_ax, 'ylim');
    if strcmp(get(h_ax, 'YScale'), 'log')
      y_rel = ( log(y) - log(YLim(1)) ) / ( log(YLim(2)) - log(YLim(1)) );
    else
      y_rel = (y-YLim(1))/(YLim(2)-YLim(1));
    end
  end

  function [x, y] = rel2absCoords(h_ax, x_rel, y_rel)
    XLim = get(h_ax, 'xlim');
    if strcmp(get(h_ax, 'XScale'), 'log')
      x = exp( x_rel * ( log(XLim(2)) - log(XLim(1)) ) + log(XLim(1)) );
    else
      x = x_rel*diff(XLim)+XLim(1);
    end
    YLim = get(h_ax, 'ylim');
    if strcmp(get(h_ax, 'YScale'), 'log')
      y = exp( y_rel * ( log(YLim(2)) - log(YLim(1)) ) + log(YLim(1)) );
    else
      y = y_rel*diff(YLim)+YLim(1);
    end
  end

  function tf = coordsWithinLimits(h_ax, x, y)
    % check if the given point (x,y) is within the limits of the axis h_ax
    XLim = get(h_ax, 'xlim');
    YLim = get(h_ax, 'ylim');
    tf = x>XLim(1) && x<XLim(2) && y>YLim(1) && y<YLim(2);
  end

  function setNewLimits(ax, xlim, ylim)
    validX = ~any(isnan(xlim)) && ~any(isinf(xlim)) && diff(xlim)>0;
    if strcmp(get(ax,'XScale'),'log')
      validX = validX && xlim(1) ~= 0;
    end
    if validX
      set(ax, 'Xlim', xlim);
    else
      if strcmp(tX.Running, 'off')
        old_color = get(ax, 'YColor');
        set(ax,'YColor','r');
        tX.TimerFcn = @(x,y)set(ax,'YColor',old_color);
        start(tX);
      end
    end
    
    validY = ~any(isnan(ylim)) && ~any(isinf(ylim)) && diff(ylim)>0;
    if strcmp(get(ax,'YScale'),'log')
      validY = validY && ylim(1) ~= 0;
    end
    if validY
      set(ax, 'Ylim', ylim);
    else
      if strcmp(tY.Running, 'off')
        old_color = get(ax, 'XColor');
        set(ax,'XColor','r');
        tY.TimerFcn = @(x,y)set(ax,'XColor',old_color);
        start(tY);
      end
    end
  end

  function i_ax = affectedAxes(currAx)
    % not speed optimized yet...
    i_ax = [];
    affAx = [];
    for i=1:numel(h_ax_linking) % iterate over linking-arrays in the cell
      if any(h_ax_linking{i} == currAx) % if current axes is part of the ith-array
        affAx = [affAx, h_ax_linking{i}]; % add axes of this array to affected axes
      end
    end
    affAx = unique(affAx);
    for i=1:numel(affAx)
      i_ax = [i_ax, find(h_ax == affAx(i))];
    end
  end

  function h_ax_linking = linkPlotyyAxes(h_ax)
    % link axes belonging to a single plotyy, all normal axes will be unlinked
    axesToBeLinked = h_ax;
    h_ax_linking = {};
    while ~isempty(axesToBeLinked)
      if isfield(getappdata(axesToBeLinked(1)), 'graphicsPlotyyPeer');
        h_PlotyyPeer = getappdata(axesToBeLinked(1), 'graphicsPlotyyPeer');
        h_ax_linking = {h_ax_linking{:}, [axesToBeLinked(1), h_PlotyyPeer]};
        axesToBeLinked = axesToBeLinked(~(axesToBeLinked==axesToBeLinked(1)) & ~(axesToBeLinked==h_PlotyyPeer));
      else
        h_ax_linking = {h_ax_linking{:}, axesToBeLinked(1)};
        axesToBeLinked = axesToBeLinked(~(axesToBeLinked==axesToBeLinked(1)));
      end
    end
  end
end