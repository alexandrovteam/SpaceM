function [vScales, vStatdists] = AutoScales( cX, cOpts )

%
% function vScales = AutoScales( cX, cOpts )
%
% Picks O(log(N)) random pts in cX, computes all distances from these points to all the other points,
% takes the average/min/max distances D, and returns a vector vScales of radii such that a ball of
% radius vScales(k) contains in average/at least/at most k*MinNperBin points.
%
% IN:
%   cX               : Dim by NumberOfPoints matrix of data points
%   [cOpts]          : structure of options with the following fields:
%                       [MinNperBin]   : a number or a # of scales by 1 vector, number of points per radial bin. 
%                                                 Default: 30.   modified by YM
%                       [StatType]     : {'Mean','Min','Max','Median'}. Take mean/min/max/median of the distances
%                                        between the random points and all points. Default: 'Mean'.
%
% OUT:
%   vScales          : column vector of scales
%   [vStatdists]     : n-th neighbor's statistical distance, added by YM
%
% USES:
%   nrsearch
%
%
% (c) Copyright 2008
% Mauro Maggioni, Duke University
%


vScales = [];

lN = size(cX,2);

if nargin<2,
    cOpts = [];
end;

if ~isfield(cOpts,'MinNperBin'), cOpts.MinNperBin=30;   end;
if ~isfield(cOpts,'StatType'),   cOpts.StatType='Mean'; end;

if length(cOpts.MinNperBin) == 1,
    if cOpts.MinNperBin>lN-1,
        cOpts.MinNperBin = lN-1;
    end;
else
    cOpts.MinNperBin = sort(cOpts.MinNperBin);
    if cOpts.MinNperBin(end) > lN
       properN = length(find(cOpts.MinNperBin <= lN));
       cOpts.MinNperBin = cOpts.MinNperBin(1:properN);
    end
end

% Pick about \log(NumberOfPoints) random points
lSeedPtIdx = randperm(lN); lSeedPtIdx = lSeedPtIdx(1:min([ceil(20*log(lN)),lN]));

% Compute the distance between the seed point and all the other points
[tmp_count,tmp_idxs,dists] = nrsearch( cX, uint32(lSeedPtIdx), lN, 0, struct('ReturnAsArrays',true,'FastNNSearcher','nn_search'));

% For debug purposes only: %figure;PlotWithStd(1:size(dists,1),dists);title('Distribution of distances from random sample to all');

switch lower(cOpts.StatType)
    case 'mean'
        statdists = mean(dists,1);
    case 'min'
        statdists = min(dists,[],1);
    case 'max'
        statdists = max(dists,[],1);
    case 'median'
        statdists = median(dists,1);
    case 'minext'
        lMin = min(dists,[],1);
        statdists = union(lMin/2,lMin);
        vScales = statdists;
    case 'maxext'
        lMax = max(dists,[],1);
        statdists = union(lMax,lMax*2);
        vScales = statdists;
    case 'ext'
        lMin = min(dists,[],1);
        lMax = max(dists,[],1);
        statdists = union(union(lMin/2,lMin),union(lMax,lMax*2));
        vScales = statdists;
    otherwise
        error(sprintf('AutoScales: unknown StatType %s',cOpts.StatType));
end;
        
% Now create bins with equal number of points
if isempty(vScales),
    if length(cOpts.MinNperBin) == 1
        vScales = statdists(cOpts.MinNperBin:cOpts.MinNperBin:lN);
    else 
        vScales = statdists(cOpts.MinNperBin);
    end
end;

if nargout == 2,
    vStatdists = statdists;
end

return;
