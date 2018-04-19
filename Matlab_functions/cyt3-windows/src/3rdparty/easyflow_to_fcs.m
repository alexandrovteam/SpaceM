function easyflow_to_fcs(data, sample_names, channels)
% Don't run the function! write in the command window:easyflow_to_fcs(data,sample_names,channels)
% I can give any name that i want to the variables (it doesn't have to be
% data, sample_names, and channels. i just have to maintain the same order.
% if i choose different names write: easyflow_to_fcs(exp1,exp2,exp3)
% the new fcs files will be saved in this matlab folder.
% Check that: - the number in data and in sample_names is the same, - that
% the number of channels is the same as in (size(data{1}))- the number of
% columns in each of the cells in the "data" variable
%
% easyflow_to_fcs(data, sample_names, channels)
%
% Convert exported EasyFlow data into FCS files.
%
% data = a cell array of matrices. each cell in the array is data from one sample.
% sample_names = the name of the samples.
% channels = the name of the channels.
%
% example:
% easyflow_to_fcs(Cyto_January, Cyto_SampleNames_January, Cyto_Channels_January)

if(size(data, 2) ~= size(sample_names, 1))
	error 'number of cells in data and sample_names is different'
end

directory_name = uigetdir(pwd,'Save FCS files')

% exporting
for sample_idx = 1:length(sample_names)
	filename = [directory_name filesep sample_names{sample_idx} '.fcs'];
	fca_writefcs(filename, data{sample_idx}, channels, channels);
end

