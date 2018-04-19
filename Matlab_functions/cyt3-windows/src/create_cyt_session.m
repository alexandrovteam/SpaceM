function create_cyt_session(X, channel_names, pathname, filename)
% X - data matrix nXm
% channel_names - OPTIONAL 1Xm length cell array of strings

if (~exist('filename', 'var') || ~exist('pathname', 'var'))
    % get filename to save to
    [filename,pathname,~] = uiputfile('*.mat','Save Session');
end

if isequal(filename,0) || isequal(pathname,0)
    return;
end

sessionData = X;

gates = cell(1, 4);
gates{1, 1} = char(filename);
gates{1, 2} = 1:size(sessionData,1);

if (exist('channel_names', 'var'))
    gates{1, 3} = channel_names;
else
    gates{1, 3} = strcat({'channel '},int2str((1:size(sessionData,2)).'))';
end
gates{1, 4} = [pathname filename]; % opt cell column to hold filename

save([pathname filename], 'sessionData', 'gates'); 


end