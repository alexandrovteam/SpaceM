function [f,g,h] = plot3k(L,varargin)
% PLOT3K  Color coded 3D scatterplot
%   [f,g,h] = plot3k(L,varargin)
% Input
%   L    An (x,y,z) vector of data points (n x 3) or
%        an (n x m) matrix with m > 3 or
%        a 3 element cell array containing x, y and z vectors or matrices.
%        For cell array input x and y may be matrices the same size as z
%        or vectors of row and column locations on a grid.  If L is an n x m
%        matrix, an xy grid is created using meshgrid(1:n,1:m).
%
%   A 3D scatterplot of the data points in L is produced using the magnitude
%   of z to color the points according to their distance from the xy plane.
%   Plot3 is called once for each group of points that map to the same color.
%   The color mapping, color range, point markers and labels can be changed
%   with optional property,value pair arguments.
%
%   varargin       'property','value' pairs that modify plot characteristics.
%      Property names must be specified as quoted strings. Property names and
%      their values are not case sensitive. For each property the default value
%      is given below.  All properties are optional. 
%
%   'ColorData'    Array of color values (default = 3rd column of L)
%      The data points are colored according to the colordata array.  The
%      number of values must equal the length of L.  A row vector, column
%      vector or matrix can be specified.
%
%   'ColorRange'   Color mapping limits  (default = [min(L(:,3) max(L,:3)])
%      A two element vector specifying the values that map to the first and
%      last colors.  This is useful for generating a series of plots with
%      identical coloring.  The colormap (but not the colorbar) is flipped
%      upside down if 'ColorRange' is given as [max min] instead of [min max]. 
%
%   'Marker'       Cell array containing point marker character and size
%      Any marker from the standard plot set '+o*.xsd^v><ph' and its size.
%      The default is {'.' 6}.
%
%   'CBLabels'     Number of colorbar labels (default = 11)
%
%   'CBFormat'     Format string for printing colorbar labels
%      Any sprintf format string with a numeric conversion code.
%      The default is '%5.3g'.
%
%   'Labels'       A cell array of five label strings (default = {})
%      These strings are used as the plot title, x axis label, y axis label,
%      z axis label and colorbar title.  To skip a string specifiy ''.
%
%   'PlotProps'    A cell array of 'property','value' pairs
%      Plotprops are applied after the plot is drawn.
%      The default value is an empty cell array.
%
% Output
%   f,g,h         Handles for the figure, axes and colorbar
%
% Example
%   figure('color','white');
%   [x,y,z] = peaks(101);
%   c = gradient(z);
%   k = hypot(x,y)<3;
%   plot3k({x(k) y(k) z(k)},...
%      'ColorData',c(k),'ColorRange',[-0.5 0.5],'Marker',{'o',2},...
%      'Labels',{'Peaks','Radius','','Intensity','Lux'},...
%      'PlotProps',{'FontSize',12});
%
% Ken Garrard, North Carolina State University, 2005 - 2011
% Based on plot3c by Uli Theune, University of Alberta
%
% 2008.11.07  Updated color range and marker arguments
% 2010.12.14  Updated input L to allow x,y vectors and z matrix
% 2011.01.04  Include axes and colorbar handles as output arguments
% 2011.10.14  Optional arguments are now property,value pairs
% 2012.03.08  The NextPlot property treatment has been updated

% -- Help

% Plot3k was called without arguments, plot an example in a new figure
% and display help text
if nargin < 1
   figure('color','white');
   [x,y,z] = peaks(101);
   c = gradient(z);
   k = hypot(x,y)<3;
   plot3k({x(k) y(k) z(k)},...
      'ColorData',c(k),'ColorRange',[-0.5 0.5],'Marker',{'o',2},...
      'Labels',{'Peaks','Radius','','Intensity','Lux'},...
      'PlotProps',{'FontSize',12});
   error(['No input arguments given\n'...
          'Please consult the help text and the example plot\n'...
          '--------\n%s'],help(mfilename));
end

%-- Parse and validate input arguments

% Convert singleton cell array to its contents
while iscell(L) && length(L) == 1, L = L{1}; end

% Check for cell array location data
if iscell(L)
   
   % cell array input must have three elements
   if length(L) ~= 3
      error('Location cell array must contain three elements');
   end

   % convert cell array elements to column vectors
   x = L{1}(:);   y = L{2}(:);   z = L{3}(:);

   % make a grid from col,row vectors
   if isequal(numel(x)*numel(y),numel(z))
      [x,y] = meshgrid(x,y);

   % equal length x, y and z vectors ?
   elseif ~isequal(length(x),length(y),length(z))
      error('Location cell array elements must be same length');
   end

   % convert to n x 3 matrix
   L = [x(:) y(:) z(:)];

% Location data may be a 2d matrix
elseif ndims(L) == 2 && size(L,2) > 3 %#ok<ISMAT>

   % use pixel values for x,y
   [x,y] = meshgrid(1:size(L,2),1:size(L,1));

   % convert to n x 3 matrix
   L = [x(:) y(:) L(:)];

% Otherwise location data must have 3 columns
elseif size(L,2) ~= 3
   error('Location vector must have 3 columns');
end

% Set up property structure with default values
p.colordata   = [];              % color by magnitude of z
p.colorrange  = [];              % color range = full range of z
p.marker      = {'.' 6};         % point marker
p.cblabels    = 11;              % number of colorbar labels
p.cbformat    = '%5.3g';         % colorbar label format string
p.labels      = {};              % no labels
p.plotprops   = {};              % no additional plot properties

% Parse property value pairs, replace defaults with values specified by caller
try
   p = structrecon(pv2struct(varargin{:}),p);
catch ERR
   error('Error parsing varargin list\n??? %s',ERR.message);
end

% Color vector must be same length as x,y,z location data or empty
p.colordata = p.colordata(:);
if ~isempty(p.colordata) && (size(L,1) ~= length(p.colordata))
   error('Location vector and color data vector must be the same length');
end

% Validate marker character
if ~iscell(p.marker), mark_ch = p.marker;  mark_sz = 6;
else
    mark_ch = p.marker{1};
    if length(p.marker)>1, mark_sz = p.marker{2};
    else                   mark_sz = 6;
    end
end
marker_str = '+o*.xsd^v><ph';
if (length(mark_ch) ~= 1) || isempty(strfind(marker_str,mark_ch))
   error('Invalid marker character, select one of ''%s''', marker_str);
end
mark_sz = max(mark_sz,0.5);

% Validate label strings cell array structure
if ~iscell(p.labels)
   error('labels argument should be a cell array');
end


% -- Setup colormap, sort data by color

% Find color limits and range
if isempty(p.colordata)                  % user specified colordata ?
   p.colordata = L(:,3);                 % no, color by magnitude of z
end
if isempty(p.colorrange)                 % user specified colorrange ?
   min_c = min(p.colordata);             % no, use colordata range
   max_c = max(p.colordata);
else                                     % user specified color range
   min_c = p.colorrange(1);
   max_c = p.colorrange(2);
end
range_c = max_c - min_c;

% Get current colormap
cmap = colormap;
if range_c < 0, cmap = flipud(cmap); end
clen = length(cmap);

% Calculate color value for each point
L(:,4) = ...
   min(max(round((p.colordata-min(min_c,max_c))*(clen-1)/abs(range_c)),1),clen);

% Sort by color value
L = sortrows(L,4);


%-- Plot data points in groups by color value

% Build index vector of color transitions (last point for each color)
dLix = [find(diff(L(:,4))>0); size(L,1)];

nxtplot = get(gca,'NextPlot');     % save 'NextPlot' property value
s = 1;                             % index of 1st point in a  color group
for k = 1:length(dLix)             % loop over each non-empty color group
   plot3(L(s:dLix(k),1), ...       % call plot3 once for each color group
         L(s:dLix(k),2), ...
         L(s:dLix(k),3), ...
         mark_ch,           ...
         'MarkerSize',mark_sz, ...
         'MarkerEdgeColor',cmap(L(s,4),:), ... % same marker color from cmap
         'MarkerFaceColor',cmap(L(s,4),:), ... % for all points in group
         'Tag','plot3k');                      % tag each group
   s = dLix(k)+1;                  % next group starts at next point
   if k == 1, hold on; end         % add next group to same axes
end
set(gca,'NextPlot',nxtplot);       % restore 'NextPlot' property


% -- Add title and axes labels

fontargs = {'FontName','Arial','FontSize',10,'FontWeight','bold'};
s = repmat({''},5-length(p.labels),1);
strlab = {p.labels{:} s{:}}; %#ok<CCAT>
title (strlab{1},fontargs{:});
xlabel(strlab{2},fontargs{:});
ylabel(strlab{3},fontargs{:});
zlabel(strlab{4},fontargs{:});

% Set plot characteristics
view([-32,32]);
grid on;
set(gca,fontargs{:});

% Apply user specified plot properties
if ~isempty(p.plotprops), set(gca,p.plotprops{:}); end


% -- Format the colorbar

h    = colorbar;
nlab = abs(p.cblabels);                  % number of labels must be positive
set(h,'YLim',[1 clen]);                  % set colorbar limits
set(h,'YTick',linspace(1,clen,nlab));    % set tick mark locations

% Create colorbar tick labels based on color vector data values
tick_vals = linspace(min_c,max_c,nlab);

% Create cell array of color bar tick label strings
labels = cell(1,nlab);
for i = 1:nlab
   labels{i} = sprintf(p.cbformat,tick_vals(i));
end

% Set tick label strings
set(h,'YTickLabel',labels,fontargs{:});
title(h,strlab{5});

if nargout > 0, f = gcf; end             % return figure handle
if nargout > 1, g = gca; end             % return axes   handle
end

% Structure reconciliation with a template
function T = structrecon(S,D)

% Check arguments, must have two structures
if ~(isstruct(S) && isstruct(D))
   error('input arguments must be structures');
end
   
T     = D;             % copy the template
fname = fields(T);     % make a list of field names

% Loop over all fields in the template, copy matching values from S
for k = 1:length(fname)
   % Process matching field names in S
   if isfield(S,fname{k})
      % Is this a substructure ?
      if isstruct(T.(fname{k})) && isstruct(S.(fname{k}))
         % Recursively process the substructure
         T.(fname{k}) = structrecon(S.(fname{k}),T.(fname{k}));
      % Not a substructure, copy field value from S
      else T.(fname{k}) = S.(fname{k});
      end
   end
end
end

% Convert argument pairs to a structure
function S = pv2struct(varargin)

% No inputs, return empty structure
if isempty(varargin), S = struct(); return; end

% Need pairs of inputs
if mod(length(varargin),2)==1
   error('number of arguments must be even');
end

% Odd elements of varargin are fields, even ones are values
% Store all field names in lower case
for k = 1:2:length(varargin)
   S.(lower(varargin{k})) = varargin{k+1};
end
end
