function [count,idxs,dists,NNInfo] = nrsearch(cX, cXq, kNN, radius, opts)

%
% function [count, idxs, dists] = nrsearch(cX, cXq, NN, radius, opts)
%
% NRSEARCH is a wrapper for performing approximate nearest neighbor
% searches using any of the following:
% 1) ANN package of Mount, et al. http://www.cs.umd.edu/~mount/ANN/
% 2) nn_search functions in the OpenTSTool box (http://www.physik3.gwdg.de/tstool/)
% 3) covertrees (to be fully integrated in V2.0 of the Diffusion Geometry toolbox)
%
% IN:
%    cX         : DxN matrix of N points in R^D comprising the data set.
%                   It may also be an L row vector of indices into cXq if cX is a subset of cXq.
%    cXq        : DxL matrix of L query points in R^D. Must be of type double.
%                   If of type 'uint32', then the query points are a subset of cX, and cXq is a row vector of length L, containing the indices of the columns of cX.
%    kNN        : number of nearest neighbors to find
%    radius     : range search with this radius (and at most NN nearest neighbors). Disregarded if 0.
%    [opts]     : structure containing the following fields:
%                   [DistRestrict]     : N array of cells, the k-th cell contains the indices of the points among which
%                                        the nearest neighbors of point k are searched for.
%                                        If it is only one cell, then it is assumed to contain the *relative* (w.r.t. current point)
%                                        indices of the points among which the nearest neighbors of point k are searched for.
%                                        Indices out of range are automatically clipped.
%                   [Tolerance]        : tolerance in range search. Default: 0 (exact).
%                   [ReturnDistSquared]: obsolete
%                   [ReturnAsArrays]   : returns <idxs> and <dists> as arrays rather than as cell arrays. This is only
%                                        possible for NN searches, not for radius searches. Default: 0.
%                   [DistInfo]         : if not empty, it will be used for distance computation. It is assumed that it contains at least
%                                        the neighbors up to the requested distance Radius (and PCARadius if provided), or kNN (if provided).
%                                        It assumes the neighbors are listed in increasing distance. It's a structure containing the following fields:
%                                        In this case cX and cXq are indices, not data points.
%                                        It may be in two formats:
%                                           idxs        : N cell, the k-th cell contains the indices of the nearest neighbors of the k-th point
%                                           dists       : N cell, the k-th cell contains the distances of the nearest neighbors of the k-th point (corresponding to idxs).
%                                        or:
%                                           count       : N vector, k-th entry being number of neighbors of k-th point stored
%                                           idxs        : N by N matrix, the k-th row contaning the indices of neighbors, in order of increasing distance
%                                           dists       : N by N nmatrix, the k-th row containing the distances, in increasing order, from the k-th point
%                   [NNInfo]            : extra structure for nearest neighbor extra information.
%                                         [CoverTree]   : If cover trees are used, then this is assumed to be a cover tree for cX.
%                                         [Atria]       : If the TSTOOL nn_prepare is used, this is the atria structure of cX.
%                   [FastNNSearcher]    : 'nn_search','ANNsearch','covertree','best'
%                                           (i.e. tries nn_search first, then ANNsearch). Default: 'best'.
%                   [SaveOpts]          : save options to a global variable called nr_search_opts. Default: false;
%                   [ReuseOpts]         : reuse global options nr_search_opts. Default: false;
%                   [Verbose]           : verbosity level 0,1,2. Default: 0.
%                   [CoverTreeOpts]     : options for covertree, if covertree-based search is used (and no NNInfo.CoverTree is provided).
%                   [distancefcn]       : distance function to use in covertree.
%
% OUT:
%    count     : 1*L vector of number of neighbors
%    idxs      : L cell array (or L times NN array if opts.ReturnAsArrays==1) of neighbor indices
%    dists     : L cell array (or L times NN array if opts.ReturnAsArrays==1) of distances of the neighbors
%    NNInfo.Atria or NNInfo.CoverTree will be populated if TSTool or covertrees (resp.) were used.
%
%
% (c) Copyright Duke University, 2008-2013
% Mauro Maggioni
%
% EXAMPLE:
%   X=[cos(linspace(0,2*pi-2/500*pi,500));sin(linspace(0,2*pi-2/500*pi,500))];
%   X=X+0.001*randn(size(X));
%   [count,idxs, dists,NNInfo] = nrsearch(X, X, 3, 0)
%   [count,idxs, dists,NNInfo] = nrsearch(X, X, 3, 0,struct('ReturnAsArrays',0));
%   [count,idxs, dists,NNInfo] = nrsearch(X, uint32((1:3)'), 3, 0,struct('ReturnAsArrays',0));

persistent ISOSX
persistent nr_search_installed_CoverTrees;
persistent nr_search_installed_TSTool;
persistent nr_search_installed_ANNsearch;
if isempty(nr_search_installed_CoverTrees), nr_search_installed_CoverTrees  = (exist('covertree')==3);  end
if isempty(nr_search_installed_TSTool),     nr_search_installed_TSTool      = (exist('nn_prepare')==3); end
if isempty(nr_search_installed_ANNsearch),  nr_search_installed_ANNsearch   = (exist('ANNsearch')==3);  end

% Setup parameters
if nargin<5,    opts = []; end
if (nargin<4) || isempty(radius),    radius = 0; end

if isfield(opts,'ReuseOpts') && (opts.ReuseOpts),
    global nr_search_opts;
    nr_search_opts.ReuseOpts = opts.ReuseOpts;
    nr_search_opts.SaveOpts = opts.SaveOpts;
    opts = nr_search_opts;
else
    if ~isfield(opts,'DistRestrict'),           opts.DistRestrict = {}; end
    if ~isfield(opts,'Tolerance'),              opts.Tolerance = 0; end
    if isfield(opts,'ReturnDistSquared'),       warning('nrsearch - ReturnDistSquared option obsolete!!'); end
    if ~isfield(opts,'ReturnAsArrays'),         opts.ReturnAsArrays = 0; end
    if ~isfield(opts,'DistInfo'),               opts.DistInfo = []; end
    if ~isfield(opts,'NNInfo'),                 opts.NNInfo = []; end
    if ~isfield(opts.NNInfo,'Atria'),           opts.NNInfo.Atria = []; end
    if ~isfield(opts.NNInfo,'CoverTree'),       opts.NNInfo.CoverTree = []; end
    if isfield(opts,'XIsTransposed'),
        if opts.XIsTransposed==true, error('GraphDiffusion:XIsTransposed option is deprecated.');
        else warning on; warning('GraphDiffusion:XIsTransposed option is deprecated.'); end
    end
    if isempty(kNN),                            kNN = 0; end
    if (~isfield(opts,'FastNNSearcher')) || (isempty(opts.FastNNSearcher)),     opts.FastNNSearcher = 'best'; end
    if ~isfield(opts,'SaveOpts'),               opts.SaveOpts = false; end
    if ~isfield(opts,'ReuseOpts'),              opts.ReuseOpts = false; end
    if ~isfield(opts,'Versbose'),               opts.Verbose = 0; end
    if ~isfield(opts,'distancefcn'),            opts.distancefcn = int32(0) ; end
    if ~isfield(opts,'CoverTreeOpts') || isempty(opts.CoverTreeOpts),          
        opts.CoverTreeOpts = struct('theta',0.5,'numlevels',int32(100),'minlevel',int32(0), ...
                                        'BLOCKSIZE',int32(1024),'distancefcn',opts.distancefcn);
    end
end

idxs = []; count = []; dists = []; NNInfo = [];

% Decide on the nearest neighbor searcher to use
lUseTSTool      = false;
lUseANNsearch   = false;
lUseCoverTrees  = false;
if isfield(opts,'DistInfo') && isempty(opts.DistInfo)
    if strcmpi(opts.FastNNSearcher,'best'),
        lUseCoverTrees  = (exist('covertree')==3);
        lUseTSTool      = (exist('nn_prepare')==3);
        lUseANNsearch   = (exist('ANNsearch')==3);
    elseif strcmpi(opts.FastNNSearcher,'nn_search'),
        lUseTSTool      = nr_search_installed_TSTool;
    elseif strcmpi(opts.FastNNSearcher,'ANNsearch'),
        lUseANNsearch   = nr_search_installed_ANNsearch;
    elseif strcmpi(opts.FastNNSearcher,'covertree'),
        lUseCoverTrees  = nr_search_installed_CoverTrees;
    end
end

if ((lUseTSTool==false) && (lUseANNsearch==false) && (lUseCoverTrees==false)),
    warning('\n ***nrsearch: could not find any fast nearest neighbor package installed. See ''help nrsearch''. Aborting nearest neighbor search.');
    return;
end

lUseCoverTrees = false;                     %%%%% FOR DEBUG PURPOSES ONLY
if lUseANNsearch && lUseTSTool
    lUseANNsearch = false;
end

% Number of points in cX
lNPts = size(cX,2);
% If no cXq, then query points are all points
if isempty(cXq), cXq = int32(1:lNPts); end
lGoWithIdxs = isa(cXq,'int32');

% Construct cover tree if needed
if lUseCoverTrees,
    if isempty(opts.NNInfo.CoverTree),
        % Construct cover tree if needed
        %NNInfo.covertree_opts   = struct('theta',0.5,'numlevels',int32(100),'minlevel',int32(0),'NTHREADS',int32(feature('numCores')),'BLOCKSIZE',int32(1024));
        NNInfo.covertree_opts   = opts.CoverTreeOpts;
        NNInfo.CoverTree        = covertree_build( cX,NNInfo.covertree_opts );
    else
        NNInfo.covertree_opts   = struct('theta',opts.NNInfo.CoverTree.theta,'numlevels',int32(opts.NNInfo.CoverTree.outparams(3)),'minlevel',int32(opts.NNInfo.CoverTree.outparams(2)),'NTHREADS',int32(feature('numCores')),'BLOCKSIZE',int32(1024));
        NNInfo.CoverTree        = opts.NNInfo.CoverTree;
    end
end

if isempty(opts.DistRestrict),                                                         % Check if a set of candidate nearest neighbors is provided
    if lUseANNsearch || lUseCoverTrees,
        % Use ANNsearch
        if lGoWithIdxs, cXq = cX(:,cXq); end                                          % If cXq is a set of indices into cX, unfortunately here we need to create it as a set of points
        if radius==0,                                                                  % Find nearest neighbors with ANNsearch
            if ~lUseCoverTrees,
                if nargout>2,
                    [idxs] = ANNsearch(cX, cXq, kNN, opts.Tolerance);
                else
                    [idxs,dists] = ANNsearch(cX, cXq, kNN, opts.Tolerance);
                end
                dists = sqrt(dists);                                                    % ANNsearch returns distances squared
            else
                [idxs,dists] = covertree_nnsearch( cX, NNInfo.CoverTree, cXq, kNN );    % Find nearest neighbors with cover trees
            end
            dists = dists';
            idxs  = idxs';
        else
            if ~lUseCoverTrees,
                [~, idxs, dists] = ANNrsearch(cX, cXq, kNN, radius^2, opts.Tolerance);
                dists = sqrt(dists);
            else
                [idxs,dists] = covertree_rangesearch(cX,NNInfo.CoverTree,cXq,radius);     % Find near points with cover trees
            end
        end
    elseif lUseTSTool,                                                                   % Use nn_search
        lX = cX';
        if ~lGoWithIdxs, lXq = cXq'; else lXq=cXq; end
        % Unfortunately there's a bug in nn_search, when there is only one point
        if size(lX,1)==1,
            lOnlyOnePoint = true;
            lX = [lX;lX];
        else
            lOnlyOnePoint = false;
        end
        if (isempty(opts.NNInfo.Atria)) || (opts.NNInfo.Atria.params(1)~=size(lX,1)),
            NNInfo.Atria = nn_prepare( lX );
        else
            NNInfo.Atria = opts.NNInfo.Atria;
        end
        if radius==0,
            % Perform nearest neighbors search
            kNN = min(size(lX,1),kNN);
            [idxs,dists] = nn_search( lX, NNInfo.Atria, double(lXq), kNN,-1,opts.Tolerance );
            idxs = uint32(idxs);
            if lOnlyOnePoint,
                idxs = ones(size(idxs));
            end
        else
            % Perform range search
            [count, neighbors] = range_search( lX, NNInfo.Atria, double(lXq), radius,-1 );
            % Sort the results
            for p = 1:size(neighbors,1),
                [neighbors{p,2},lSortedIdxs] = sort(neighbors{p,2});
                neighbors{p,1} = neighbors{p,1}(lSortedIdxs);
            end
            % Convert to standard output format
            if ~lOnlyOnePoint,
                if kNN==0,
                    for k = length(count):-1:1,
                        idxs{k} = uint32(neighbors{k,1});
                        dists{k} = neighbors{k,2};
                    end
                else
                    % Truncate to nearest NN points
                    for k = length(count):-1:1,
                        count(k) = min(count(k),kNN);
                        idxs{k} = uint32(neighbors{k,1}(1:count(k)));
                        dists{k} = neighbors{k,2}(1:count(k));
                    end
                end
            else
                for k = length(count):-1:1,
                    idxs{k} = uint32(1);
                    dists{k} = neighbors{k,2};
                end
            end
            clear neighbors;
        end
    end
else
    % A list of candidate nearest neighbors was provided
    for k = lNPts:-1:1,
        if (iscell(opts.DistRestrict)) && (length(opts.DistRestrict)==lNPts),
            % A NN list is provided for each point
            lCurIdxs = opts.DistRestrict{k};
        else
            % Only one NN list, interpret as relative to the current point
            lCurIdxs = k+opts.DistRestrict;
            % Clip indices out of range
            lCurIdxs((lCurIdxs<1) | (lCurIdxs>lNPts)) = [];
        end
        % Find the neighbors
        if radius==0,
            if ~lUseTSTool,
                [idxs_tmp,dists_tmp] = ANNsearch( cX(:,lCurIdxs), cX(:,k), kNN, opts.Tolerance );
                dists_tmp = sqrt(dists_tmp);
            else
                lX = cX(:,lCurIdxs)';
                if isempty(opts.NNInfo.Atria),          % MM: Not correct if Atira is not for the restricted indices
                    NNInfo.Atria = nn_prepare( lX );
                else
                    NNInfo.Atria = opts.NNInfo.Atria;
                end
                [idxs_tmp,dists_tmp] = nn_search( lX, NNInfo.Atria, cX(:,k), kNN, opts.Tolerance );
                if opts.ReturnDistSquared,
                    for k = 1:size(idxs_tmp,1);dists_tmp(k,:) = dists_tmp(k,:).^2;end
                end
            end
        else
            if ~lUseTSTool,
                [count_tmp,idxs_tmp,dists_tmp] = ANNrsearch( cX(:,lCurIdxs), cX(:,k), kNN, radius^2, opts.Tolerance );
                dists_tmp = sqrt(dists_tmp);
            else
                lX = cX(:,lCurIdxs)';
                if isempty(opts.NNInfo.Atria) || (opts.NNInfo.Atria.params(1)~=size(lX,1)),
                    NNInfo.Atria = nn_prepare( lX );
                else
                    NNInfo.Atria = opts.NNInfo.Atria;
                end
                [count_tmp, neighbors_tmp] = range_search( lX, NNInfo.Atria, cX(k,:)', radius );
                for k = length(count_tmp):-1:1,
                    idxs_tmp{k}  = neighbors_tmp{k,1};
                    dists_tmp{k} = neighbors_tmp{k,2};
                end
                clear neighbors_tmp;
            end
        end
    end
    if ~opts.ReturnAsArrays,
        idxs{k} = uint32(lCurIdxs(idxs_tmp));
        dists{k} = dists_tmp;
    else
        idxs(k,:) = uint32(lCurIdxs(idxs_tmp));
        if argout>2,
            dists(k,:) = dists_tmp;
        end
    end
end


% Format outputs as arrays or cell, as desired
if ~opts.ReturnAsArrays,
    if ~iscell(idxs),
        idxs  = mat2cell(idxs,ones(size(idxs,1),1),size(idxs,2));
        dists  = mat2cell(dists,ones(size(dists,1),1),size(dists,2));
    end
else
    if iscell(idxs),
        idxs_2 (lNPts,lNPts) = uint32(0);
        dists_2(lNPts,lNPts) = single(Inf);
        for k = 1:length(idxs),
            idxs_2(k,1:length(idxs{k}))   = uint32(idxs{k});
            dists_2(k,1:length(dists{k})) = dists{k};
        end
        idxs = uint32(idxs_2); dists = dists_2; clear idxs_2 dists_2;
    end
end


% Build count
if (~iscell(idxs)) && (size(idxs,1)==1),
    count = ones(1,length(idxs),'uint32');
else
    if ~opts.ReturnAsArrays,
        for k = length(idxs):-1:1,
            count(k) = length(idxs{k});
        end
    else
        count = size(idxs,2)*ones(1,size(idxs,1));
    end
end

if opts.SaveOpts,
    global nr_search_opts
    if exist('NNInfo'),
        opts.NNInfo = NNInfo;
    end
    nr_search_opts = opts;
end


return;