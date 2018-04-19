function Startup_DiffusionGeometry(BaseDirectory)

% Startup_DiffusionGeometry adds the various directories used by the diffusion geometry code
% to the current path.  Startup_DiffusionGeometry will assume the current directory is the
% base directory unless another directory is specified


fprintf('Startup_DiffusionGeometry.m: setting diffusion geometry paths ... \n');

if nargin==0
    Prefix  = [pwd filesep];
else
    Prefix  = [BaseDirectory filesep];
end;

% choose your nearest neighbors code by changing this line
appendpath(([Prefix]));
appendpath(([Prefix 'Diffusions']));
appendpath(([Prefix 'NearestNeighbors']));
appendpath(([Prefix 'Utils']));
appendpath(([Prefix 'Examples']));
appendpath(([Prefix 'CoverTrees']));

fprintf('Startup_DiffusionGeometry.m: disabling case sensitivity warning ... \n');
warning('off','MATLAB:dispatcher:InexactMatch');

function appendpath(string)

fprintf('\t%s\\ \n', string);
addpath(genpath(string));

return;

