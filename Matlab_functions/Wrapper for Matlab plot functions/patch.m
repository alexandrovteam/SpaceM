function h = patch(varargin)
%% wrapper for the original patch function
currentFolder = pwd; % save current folder
cd([matlabroot '\toolbox\matlab\specgraph']) %go to matlab folder
try
  h = patch(varargin{:}); % call original function
  akZoom();
  cd(currentFolder) % go back to current folder
catch err
  cd(currentFolder) % go back to current folder
  rethrow(err)
end

% suppress output if not needed
if nargout == 0
  clear h;
end