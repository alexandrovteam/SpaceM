function [n_bColors, n_sColors] = FCmap(bmap, smap, bColors, sColors, varargin)
% bColors and sColors is a representation of fold-change in channel smoothed over map
% bMap is an nxd basal map.
% sMap is an mxd stim map. 
% bChannel and sChannel are the expression level vector respectively.  
% bColors and sColors are the smooth differences in the means of expression
% of the KNN for each point

% set defaults 
k = 8;
thresh = 0;
stim = 0;

i = 1;
while i<=length(varargin)
    if strcmp(varargin{i},'K')
        k=varargin{i+1}; 
    elseif strcmp(varargin{i},'thresh')
        thresh=varargin{i+1}; i=i+1;
    elseif strcmp(varargin{i},'stim')
        stim=varargin{i+1}; i=i+1;
    end
    i=i+1;
end


%pre-process: cap high outliers in channel
% b_idx = find_outliers_Thompson(bColors);
% s_idx = find_outliers_Thompson(sColors);
% clean_b_channel = bColors;
% clean_b_channel(b_idx) = [];
% clean_s_channel = sColors;
% clean_s_channel(s_idx) = [];
% bColors(bColors<min(clean_b_channel)) = min(clean_b_channel);
% sColors(sColors<min(clean_s_channel)) = min(clean_s_channel);
% bColors(bColors>max(clean_b_channel)) = max(clean_b_channel);
% sColors(sColors>max(clean_s_channel)) = max(clean_s_channel);

n = size(bmap,1);

hwaitbar = waitbar(0,'Build KNN graph from basal to basal');
bbIDX = knnsearch_fast( bmap, bmap, k);
waitbar(1/4,hwaitbar, 'Build KNN graph from basal to stim');
bsIDX = knnsearch_fast( bmap, smap, k);

n_bColors = mean( sColors( bsIDX(:, :) ), 2 ) - mean( bColors( bbIDX(:, :) ), 2 );
% sColors = 0;

if stim
    waitbar(2/4,hwaitbar, 'Build KNN graph from stim to basal');
    sbIDX = knnsearch_fast( smap, bmap, k);
    waitbar(3/4,hwaitbar, 'Build KNN graph from stim to stim');
    ssIDX = knnsearch_fast( smap, smap, k);
    waitbar(4/4,hwaitbar, 'Compute difference of means by k-community');

    n_sColors = mean( sColors( ssIDX(:, :) ), 2 ) - mean( bColors( sbIDX(:, 2:end) ), 2 );
else
    n_sColors = sColors*0;
end


close(hwaitbar);
end