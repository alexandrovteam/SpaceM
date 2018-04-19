% get current path
[curr_path, ~, ~] = fileparts(mfilename('fullpath'));
curr_path = [curr_path filesep];

% get all subdirectories
path_with_subdirectories = genpath(curr_path );
addpath( path_with_subdirectories );
savepath;
cyt


