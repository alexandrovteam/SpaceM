function [AX,H1,H2] = plotyy(varargin)
%% wrapper for the original plotyy function
%% avoid name conflict warning with: "warning off MATLAB:dispatcher:nameConflict"
currentFolder = pwd; % save current folder
cd([matlabroot '\toolbox\matlab\graph2d']) %go to matlab folder
try
  [AX,H1,H2] = plotyy(varargin{:}); % call original function
  plotyy_fixZoom(AX); % this function fixes problems with plotyy tickmarks (not related to akZoom)
  akZoom();
  cd(currentFolder) % go back to current folder
catch err
  cd(currentFolder) % go back to current folder
  rethrow(err)
end

% suppress output if not needed
if nargout == 0
  clear AX;
end