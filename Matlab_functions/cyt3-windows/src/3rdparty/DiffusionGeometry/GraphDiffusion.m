function G = GraphDiffusion(Points, Radius, Options)

%
% function G = GraphDiffusion(Points, Radius, Options)
%
% GraphDiffusion builds a diffusion operator on a cloud of points embedded in
% a metric space.  A matrix of graph weights for the point cloud is
% constructed by building a local bump function centered at each point and
% then the weight matrix is normalized to form the diffusion operator T.
%
% GraphDiffusion builds the weight matrix by constructing each column as a local
% bump function.  For each point x, the points within Radius are found. The
% weight
%
%    exp( - d(x,y) * d(y,x) / (epsilon_x * epsilon_y) )
%
% is then assigned to the pair (x,y).  Note that because of autotuning,
% d(x,y) is not necessarily equal to d(y,x) ... thus the product.  Also
% note that the constant epsilon can vary from point to point.
%
% The matrix of weights is then normalized, in one of a number of possible
% ways, to form a final "diffusion" operator T.  Normally, either a symmetric
% Markov or symmetric Laplace-Beltrami normalization is used in order to make
% eigensolving
%
%
% IN:
%
%    Points                 : DxN matrix specifying N points in R^D
%                             If it empty, DistInfo must be properly given.
%    Radius                 : radius of range search. If Radius==0, uses kNN only.
%                             It can be a vector of length K, in which case K graphs will be constructed, one for each value of Radius.
%    Options                : structure for specifying algorithm options which can contain
%                               zero or more of the following fields:
%
%      [kEigenVecs]         : number of eigenvectors of Laplacian to compute. Default: 0.
%      [Display]            : when this is zero the only output will be error messages
%      [Epsilon]            : either a scalar or a vector of length K (same as Radius) giving time-step constants
%                               as in the formula above. If 0, uses binary weights on edges.
%      [kNN]                : truncate the neighbors to the k closest points.
%      [kNNAutotune]        : rescales d(x,y) based on the distance to the kNNAutotune nearest point: divides d(x,y) by the distance
%                             between x and the kNNAutotune-th nearest neighbor of x.
%      [NNMaxDim]           : dimension of the space to project onto for nearest neighbor computation. Default: 0 (no projection).
%      Normalization        : method for normalizing the matrix of weights
%         'bimarkov'            force row and column sums to be 1
%         'markov'              force row sums to be 1
%         'smarkov'             symmetric conjugate to markov
%         'beltrami'            Laplace-Beltrami normalization ala Coifman-Lafon
%         'sbeltrami'           symmetric conjugate to beltrami
%         'FokkerPlanck'        Fokker-Planck normalization
%         'sFokkerPlanck'       symmetric conjugate to Fokker-Planck normalization
%      [Distance]           : distance to use
%         'Euclidean'       : default
%         Otherwise, a combination of the following strings (case sensitive):
%           'PCA'                 : project onto local tangent plane
%           'direct' or 'inverse' : local pca or its inverse
%           'cov'                 : normalize by number of local points
%           'noS'                 : no scaling by the singular values
%      [PCARadius]          : radius for the computation of local PCA. it can be a vector (matching Radius). Default: Radius.
%      [PCAtruncD]          : only use the first D principal axes for computing local PCA distance. Default: D is maximum possible.
%      [ExtraEdges]         : an m by 2 matrix of m edges to which extra weight should be added.
%      [WeightFcn]        : handle to a function that convert distances to weights (similarities). Default: uses exp. weights
%                           as described above.
%      [DistInfo]         : if not empty, it will be used for distance computation. It is assumed that it contains at least
%                         the neighbors up to the requested distance Radius (and PCARadius if provided), or kNN (if provided).
%                         It assumes the neighbors are listed in increasing distance. It's a structure containing the following fields:
%           count       : N vector, the k-th entry is the number of neighbors of the n-th point
%           idxs        : N cell, the k-th cell contains the indices of the nearest neighbors of the k-th point
%           [dists]     : N cell, the k-th cell contains the distances of the nearest neighbors of the k-th point (corresponding to idxs).
%                         Should be sorted in incrasing magnitude.
%           [CoverTree] : covertree on the N points.
%      [DistRestrict]     : N array of cells, the k-th cell contains the indices of the points among which
%                         the nearest neighbors of point k are searched for.
%      [Symmetrization]   : how to symmetrize the weight matrix:
%                         'W.Wt'     : replaces W by W.*W'
%                         'WWt'      : replaces W by W*W'
%                         'WtW'      : replaces W by W'*W
%                         'sqrtW.Wt'     : replaces W by sqrt(W.*W' (entrywise sqrt)
%                         'sqrtWWt'      : replaces W by sqrt(W*W') (entrywise sqrt)
%                         'sqrtWtW'      : replaces W by sqrt(W'*W) (entrywise sqrt)
%                         'W+Wt'     : replaces W by W+Wt. This is the default value
%                         'none'     : no symmetrization
%      [DontReturnDistInfo]: memory-saving option: DistInfo is not returned in G. Default: 0 (distinfo is returned).
%      [NoSelfLoop]        : Enforcing W(i,i)=0 to remove self-loops. Default: 0.
%      [BinaryWeights]     : If Epsilon = 0, a symmetrization may distort binary weights on edges.
%                           This option recovers binary weights on edges. Default: 0.
%      [ReturnAdjustedDistInfo]: If KNNAutotune or PCA is used, it changes the original distances. To watch the change,
%                              it reports as adjusteddists. Default: 0.
%      [W]                 : matrix of weights. If provided, only performs normalizations and eigenvector computation. Default: [].
%      [FastNNSearcher]    : passed to nrsearch for choosing the fast nearest neighbor searcher. Default: [];
%      [CoverTreeOpts]     : passed to nrsearch, options for covertree construction. See there for details.
%
%    The default options are:
%
%       Options.Display = 1;
%       Options.Epsilon = 1.0;
%       Options.Normalization = 'beltrami';
%
% OUT:
%
%    G      : Graph structure containing the following fields:
%               T         : NxN sparse matrix giving the normalized diffusion operator
%               W         : NxN sparse matrix of weights
%               P         : in the case of symmetric normalizations, this is the NxN diagonal
%                            matrix which relates the nonsymmetric version to the symmetric
%                           form via conjugation; symmetric = P^(-1)*nonsymetric*P
%               DistInfo  : structure containing distance info:
%                           count    :   N vector of counts of points near each point
%                           idxs     :   N cell of indices of the nearest points to each point
%                           dists    :   N cell of distances between the nearst points and each point
%                           [adjusteddists]: N cell of distances adjusted by KNNAutotune or PCA.
%                           [CoverTree] : covertree data structure if such was used to perform range searches.
%               Autotune  : returned only if kNNAutotune specified, the i-th entry is the distance between the i-th point and its kNNAutotune-th neighbor
%               Opts      : structure of options used to construct the graph.
%               EigenVecs : kEigenVecs eigenvectors of T, if required, one per column.
%               EigenVals : kEigenVecs eigenvalues of T, if required, column vector.
%
%
%
% Example:
%   X=[cos(linspace(0,2*pi-2/500*pi,500));sin(linspace(0,2*pi-2/500*pi,500))];
%   G = GraphDiffusion(X, 0.05, struct('Normalization','bimarkov'));

%
%
% (c) Copyright Yale University, 2006, James C Bremer Jr.
% (c) Copyright Duke University, 2007, Mauro Maggioni
%

% initialize output variables in case of error return
T = [];     W = [];     D = [];

%
% process options
%
G.Opts.kEigenVecs = 0;
G.Opts.Display = 1;
G.Opts.Epsilon = 1.0*ones(1,length(Radius));
G.Opts.Normalization = 'beltrami';
G.Opts.Distance = 'Euclidean';
G.Opts.WeightFcn = [];
G.Opts.PCARadius = Radius;
G.Opts.PCAtruncD = Inf;
G.DistInfo = [];
G.Opts.kNN = Inf;
G.Opts.kNNAutotune = 0;
G.Opts.NNMaxDim = 0;
G.Opts.Symmetrization = 'W+Wt';
G.Opts.Radius = Radius;
G.Opts.DistRestrict = {};
G.Opts.DontReturnDistInfo = 0;
G.Opts.NoSelfLoop=0;
G.Opts.BinaryWeights=0;
G.Opts.ReturnAdjustedDistInfo=0;
G.Opts.ExtraEdges = [];
G.Opts.FastNNSearcher = [];
G.Opts.W = [];
G.Opts.CoverTreeOpts = [];
if isfield(Options,'XIsTransposed'),
    if Options.XIsTransposed==true, error('GraphDiffusion:XIsTransposed option is deprecated.');
    else warning on; warning('GraphDiffusion:XIsTransposed option is deprecated.'); end
end

DistInfo = [];
lDistInfoProvided = false;

if nargin<3,    Options = struct(); end

% Go through the options specified
fn = fieldnames(Options);
for j=1:length(fn)
    name = fn{j};
    value = getfield(Options, name);
    
    if     strcmpi(name,'Display')          G.Opts.Display = value;
    elseif strcmpi(name,'kEigenVecs')       G.Opts.kEigenVecs = value;
    elseif strcmpi(name,'Epsilon')          G.Opts.Epsilon = value;
    elseif strcmpi(name,'Normalization')    G.Opts.Normalization = value;
    elseif strcmpi(name,'Distance')         G.Opts.Distance = value;
    elseif strcmpi(name,'PCARadius')        G.Opts.PCARadius = value;
    elseif strcmpi(name,'PCAtruncD')        G.Opts.PCAtruncD = value;
    elseif strcmpi(name,'DistInfo')         lDistInfoProvided = true; DistInfo = value;
    elseif strcmpi(name,'kNN')              G.Opts.kNN = value;
    elseif strcmpi(name,'kNNAutotune')      G.Opts.kNNAutotune = value;
    elseif strcmpi(name,'NNMaxDim')         G.Opts.NNMaxDim = value;
    elseif strcmpi(name,'DistRestrict')     G.Opts.DistRestrict = value;
    elseif strcmpi(name,'Symmetrization')   G.Opts.Symmetrization = value;
    elseif strcmpi(name,'DontReturnDistInfo') G.Opts.DontReturnDistInfo = value;
    elseif strcmpi(name,'WeightFcn')        G.Opts.WeightFcn = value;
    elseif strcmpi(name,'NoSelfLoop')       G.Opts.NoSelfLoop = value;
    elseif strcmpi(name,'BinaryWeights')    G.Opts.BinaryWeights = value;
    elseif strcmpi(name,'ReturnAdjustedDistInfo')  G.Opts.ReturnAdjustedDistInfo=value;
    elseif strcmpi(name,'ExtraEdges')       G.Opts.ExtraEdges = value;
    elseif strcmpi(name,'W')                G.Opts.W = value;
    elseif strcmpi(name,'FastNNSearcher')   G.Opts.FastNNSearcher = value;
    elseif strcmpi(name,'CoverTreeOpts')    G.Opts.CoverTreeOpts = value;
    else   fprintf('GraphDiffusion.m: invalid option "%s" ignored.\n', name);
    end
end

if ~isempty(Points)
    [Dim,N] = size(Points);
elseif (lDistInfoProvided) && (isfield(DistInfo,'idxs')),
    N=length(DistInfo.idxs);
elseif isempty(G.Opts.W),
    error('Points or DistInfo must be given.');
end

if G.Opts.Display > 1
    fprintf('Options: \n');
    fprintf('\tDisplay = %d\n', G.Opts.Display);
    fprintf('\tEpsilon = %g\n', G.Opts.Epsilon);
    fprintf('\tNormalization = %s\n', G.Opts.Normalization);
end


if isempty(G.Opts.W),
    %
    % Perform range/nearest neighbor searches on the input set
    %
    if G.Opts.Display; fprintf('performing range searches ... '); TIME = clock; end
    
    % Decide whether to do radius searches or nn searches
    G.Opts.DoRadiusSearch = (length(G.Opts.Radius)>1) | (G.Opts.Radius>0);
    if G.Opts.DoRadiusSearch,
        % Correct Epsilon to have the same length as Radius.
        if length(G.Opts.Epsilon)~=length(Radius),
            G.Opts.Epsilon = G.Opts.Epsilon(1)*ones(1,length(Radius));
            fprintf('\nGraphDiffusion:Warning:Value of Epsilon changed to match length of Radius parameters.');
        end
    else
        % Correct Epsilon to have the same length as kNN.
        if length(G.Opts.Epsilon)~=length(G.Opts.kNN),
            G.Opts.Epsilon = G.Opts.Epsilon(1)*ones(1,length(G.Opts.kNN));
            fprintf('\nGraphDiffusion:Warning:Value of Epsilon changed to match length of kNN parameters.');
        end
    end
    
    if isempty(DistInfo) || ~isfield(DistInfo,'idxs'),
        lOptsTemp.kNN = G.Opts.kNN;
        lOptsTemp.DistRestrict = G.Opts.DistRestrict;
        lOptsTemp.FastNNSearcher = G.Opts.FastNNSearcher;
        lOptsTemp.CoverTreeOpts = G.Opts.CoverTreeOpts;
        
        if isfield(DistInfo,'CoverTree'),
            lOptsTemp.NNInfo.CoverTree = DistInfo.CoverTree;
        else
            if (G.Opts.NNMaxDim>0) && (G.Opts.NNMaxDim<Dim) && (G.Opts.NNMaxDim<N),
                lPointsProj = Points;
                lPointsProj = bsxfun(@minus,Points,mean(lPointsProj,2));
                [U,~,~] = randPCA(lPointsProj,G.Opts.NNMaxDim);
                lPointsProj = U'*lPointsProj;
                PointsOld = Points;
                Points = lPointsProj;
            end
        end
        
        if G.Opts.DoRadiusSearch,       % Do radius search
            [DistInfo.count, DistInfo.idxs, DistInfo.dists,DistInfo.NNInfo] = nrsearch(Points, [], [], max([G.Opts.Radius,G.Opts.PCARadius]), lOptsTemp);
        else                            % Do nearest neighbors search
            [DistInfo.count, DistInfo.idxs, DistInfo.dists,DistInfo.NNInfo] = nrsearch(Points, [], max(G.Opts.kNN), [], lOptsTemp);
            DistInfo.count = max(G.Opts.kNN);
        end
        if (G.Opts.NNMaxDim>0) && (G.Opts.NNMaxDim<Dim) && (G.Opts.NNMaxDim<N) && exist('PointsOld')
            Points = PointsOld;
            clear PointsOld lPointsProj;
        end
    end
    
    if G.Opts.Display; fprintf('%g seconds\n', etime(clock, TIME)); end
    
    %
    % Compute the weights and associated operators
    %
    if G.Opts.Display; fprintf('computing weights...'); TIME = clock; end
    
    % Number of graphs to be built
    if G.Opts.DoRadiusSearch,    lNumberOfGraphs = length(G.Opts.Radius);
    else                         lNumberOfGraphs = length(G.Opts.kNN);
    end
else
    if iscell(G.Opts.W),
        lNumberOfGraphs = length(G.Opts.W);
    else
        lNumberOfGraphs = 1;
    end
    
    N = size(G.Opts.W,1);
end

% Construct the graphs and operators
for k = 1:lNumberOfGraphs,
    if isempty(G.Opts.W),
        % Find the neighbors of each point
        for i = N:-1:1,
            if G.Opts.DoRadiusSearch,
                % If there is more than one radius, or distinfo was provided externally (so we do not with what radius it was computed),
                % we need to restrict the radius to Radius(k)
                if iscell(DistInfo.dists),    lIdxs = find(DistInfo.dists{i}<=G.Opts.Radius(k));
                else                          lIdxs = find(DistInfo.dists(i,:)<=G.Opts.Radius(k));
                end
                % We restrict to the kNN nearest neighbors
                if length(lIdxs)>G.Opts.kNN,  lIdxs = lIdxs(1:G.Opts.kNN);
                end
                % Save the local distance information for point i
                if iscell(DistInfo.dists),
                    dists{i} = DistInfo.dists{i}(lIdxs);
                    idxs{i}  = DistInfo.idxs {i}(lIdxs);
                else
                    dists{i} = DistInfo.dists(i,lIdxs);
                    idxs{i}  = DistInfo.idxs (i,lIdxs);
                end
                count(i) = length(lIdxs);
            else
                % If NN's search
                if iscell(DistInfo.dists),
                    if (length(G.Opts.kNN)>1) | (lDistInfoProvided==true),
                        try
                            lIdxs = 1:min([length(DistInfo.idxs{i}),G.Opts.kNN(k)]);
                        catch
                            1,
                        end
                        idxs{i}  = DistInfo.idxs{i}(lIdxs);
                        dists{i} = DistInfo.dists{i}(lIdxs);
                    else
                        idxs{i}  = DistInfo.idxs{i};
                        dists{i} = DistInfo.dists{i};
                    end
                else
                    if (length(G.Opts.kNN)>1) | (lDistInfoProvided==true),
                        lIdxs = 1:min([DistInfo.count(i),G.Opts.kNN(k)]);
                        idxs{i}  = DistInfo.idxs(i,lIdxs);
                        dists{i} = DistInfo.dists(i,lIdxs);
                    else
                        idxs{i}  = DistInfo.idxs(i,:);
                        dists{i} = DistInfo.dists(i,:);
                    end
                end
                count(i) = length(idxs{i});
            end
        end
        % Do the autotuning of distances if requested
        if G.Opts.kNNAutotune > 0
            if G.Opts.Display; fprintf('autotuning distances ... '); TIME = clock; end
            for j=N:-1:1,
                temp = sort(dists{j}, 'ascend');
                lMaxTempIdxs = min([G.Opts.kNNAutotune,length(temp)]);
                if (lMaxTempIdxs==0) || (temp(lMaxTempIdxs)==0),        dists{j} = 0;
                else                                                    dists{j} = dists{j}./temp(lMaxTempIdxs);
                end
                G.Autotune(j) = single(temp(lMaxTempIdxs));
            end
            if G.Opts.Display; fprintf('%g seconds\n', etime(clock, TIME)); end
        end
        % Do the local PCA warping of the distances
        if ~isempty(strfind(G.Opts.Distance,'PCA')),
            if G.Opts.Display; fprintf('adjusting by local PCA ... '); TIME = clock; end
            % Compute the local PCA matrices
            for i = N:-1:1,
                % Find the nearest neighbors to be used for the local PCA computation
                if G.Opts.DoRadiusSearch,
                    % MM: TODO, when DistInfo contains matrices
                    pca_idxs  = DistInfo.idxs{i}(find(DistInfo.dists{i}<=G.Opts.PCARadius(k)));
                    pca_count = length(pca_idxs);
                else
                    % Use k nearest neighbors to compute local pca
                    % MM: TODO, when DistInfo contains matrices
                    pca_idxs  = DistInfo.idxs{i};
                    if iscell(DistInfo.count),    pca_count = DistInfo.count{i};
                    else                          pca_count = DistInfo.count;
                    end
                end
                % Center the points
                lCenteredPcaPts = bsxfun(@minus,Points(:,pca_idxs),-Points(:,i));
                % Compute local PCA
                [lU{i},lS{i}] = svd(lCenteredPcaPts,0);lS{i}=diag(lS{i});
                % Fix case in which there were less points than the dimension
                if size(lCenteredPcaPts,2)<size(lCenteredPcaPts,1),         % If there are less points than dimensions
                    lU{i} = [lU{i},zeros(size(lU{i},1),size(lCenteredPcaPts,1)-size(lCenteredPcaPts,2))];
                    lS{i} = [lS{i};zeros(size(lCenteredPcaPts,1)-size(lCenteredPcaPts,2),1)];
                end
                % Truncate the PCAs if requested
                if G.Opts.PCAtruncD<length(lS{i}),
                    lS{i}(G.Opts.PCAtruncD+1:length(lS{i}))=0;
                end
                % Compute distances between points based on the local PCA
                lCenteredRadiusPts = bsxfun(@minus,Points(:,idxs{i}),Points(:,i));
                if ~isempty(strfind(G.Opts.Distance,'cov')),
                    try
                        lS{i} = lS{i}/double(pca_count);
                    catch
                        fpritnf('Ops!');
                    end
                end
                if ~isempty(strfind(G.Opts.Distance,'noS')),
                    dists{i} = sqrt(sum((lU{i}*lCenteredRadiusPts).^2,1)');
                elseif ~isempty(strfind(G.Opts.Distance,'direct')),
                    dists{i} = sqrt(sum((diag(lS{i})*lU{i}*lCenteredRadiusPts).^2,1)');
                elseif ~isempty(strfind(G.Opts.Distance,'inverse')),
                    try
                        dists{i} = sqrt(sum((diag(lS{i}.^(-1))*lU{i}*lCenteredRadiusPts).^2,1)');
                    catch
                        fprintf('Something wrong.');
                    end
                end
            end
            if G.Opts.Display; fprintf('%g seconds\n', etime(clock, TIME)); end
        end
        
        if G.Opts.ReturnAdjustedDistInfo && (G.Opts.kNNAutotune || ~isempty(strfind(G.Opts.Distance,'PCA')))
            DistInfo.adjusteddists=dists';
        end
        
        % Build the weight matrix
        if G.Opts.Display; fprintf('forming weight matrix ... '); TIME = clock; end
        % First form the distance matrix
        rowidxs = zeros(sum(count), 1,'double');
        colidxs = zeros(sum(count), 1,'double');
        distsmat = zeros(sum(count), 1);
        index = 1;
        location = 1;
        for j=1:length(count),
            rowidxs(location:(location+count(j)-1))     = idxs{index};
            colidxs(location:(location+count(j)-1))     = index;
            distsmat(location:(location+count(j)-1))    = dists{index};
            location                                    = location+count(j);
            index                                       = index+1;
        end
        
        if G.Opts.Epsilon(k)>0,     W = sparse(rowidxs, colidxs,  distsmat, N, N, sum(count));
        else                        W = sparse(rowidxs, colidxs,  ones(size(distsmat)), N, N, sum(count));
        end
        clear idxs rowidxs colidxs dists count;
        
        if ~isempty(G.Opts.ExtraEdges),
            W2 = sparse(G.Opts.ExtraEdges(:,1),G.Opts.ExtraEdges(:,2), 1, N,N, size(G.Opts.ExtraEdges,1));
            W  = (W+W2)/2;
        end
        
        % Now take d(x,y) = d(x,y)*d(y,x) in order to symmetrize
        if ~isempty(strfind(G.Opts.Symmetrization,'W.Wt')),        W = W .* W';
        elseif ~isempty(strfind(G.Opts.Symmetrization,'WWt')),     W = W*W';
        elseif ~isempty(strfind(G.Opts.Symmetrization,'WtW')),     W = W'*W;
        elseif ~isempty(strfind(G.Opts.Symmetrization,'W+Wt')),    W = W + W';
        end
        if ~isempty(strfind(G.Opts.Symmetrization,'sqrt')),        W = sqrt(W);
        end
        
        if isempty(G.Opts.WeightFcn),
            % Adjust with epsilon
            if G.Opts.Epsilon(k)>0,
                [i j s] = find(W);
                nnzW=nnz(W);
                clear W;
                if G.Opts.NoSelfLoop,
                    s = s./(G.Opts.Epsilon(k)^2);
                    W = sparse(i, j, exp(-s), N, N, nnzW);
                else
                    i = [i; (1:N)'];
                    j = [j; (1:N)'];
                    s = [s./(G.Opts.Epsilon(k)^2); zeros(N,1)];
                    W = sparse(i, j, exp(-s), N, N, nnzW+N);
                end
            else
                if G.Opts.NoSelfLoop,
                    W = W - sparse(diag(diag(W)));
                end
                if G.Opts.BinaryWeights,
                    [i, j]=find(W);
                    nnzW=nnz(W);
                    W=sparse(i, j, 1, N, N, nnzW);
                end
            end
        else
            for k = 1:N,
                lIdxs = setdiff(find(W(k,:)),k);
                if ~isempty(lIdxs),
                    lW = feval( G.Opts.WeightFcn, Points, [k,lIdxs], [] );
                    W(k,lIdxs) = lW; W(k,k) = 1;
                end
            end
        end
        
        if G.Opts.Display; fprintf('%g seconds\n', etime(clock, TIME)); end
    else
        W = G.Opts.W;
    end
    P = [];
    
    %
    % Form the diffusion operator
    %
    if ~strcmpi(G.Opts.Normalization, 'bimarkov')
        if G.Opts.Display; fprintf('normalizing to form operator '); TIME = clock; end
    end
    
    % Create D
    DInv = sum(W,2);
    lNonZeroIdxs = find(DInv~=0);
    DInv(lNonZeroIdxs) = DInv(lNonZeroIdxs).^(-1);
    
    % Go through the various normalizations
    if strcmpi(G.Opts.Normalization, 'bimarkov')
        T = Bimarkov(W);
    elseif strcmpi(G.Opts.Normalization, 'smarkov')
        if G.Opts.Display; fprintf('(symmetric markov) ... '); end
        
        D = sparse(1:N, 1:N, sqrt(DInv), N, N, N);
        P = D;
        T = D*W*D;
        
        T = (T'+T)/2;  % iron out numerical wrinkles
        
    elseif strcmpi(G.Opts.Normalization, 'markov')
        if G.Opts.Display; fprintf('(markov) ... '); end
        T = sparse(1:N, 1:N, DInv, N, N, N) * W;
        
    elseif strcmpi(G.Opts.Normalization, 'sbeltrami')
        
        if G.Opts.Display;fprintf('(symmetric beltrami) ... '); end
        
        % beltrami normalization
        P = sparse(1:N, 1:N, DInv, N, N, N);
        K = P*W*P;
        
        DInv = sum(K,2);
        lNonZeroIdxs = find(DInv~=0);
        DInv(lNonZeroIdxs) = sqrt(DInv(lNonZeroIdxs).^(-1));
        
        D = sparse(1:N, 1:N, DInv, N, N, N);
        P = D;
        T = D*K*D;
        
        T = (T'+T)/2;  % iron out numerical wrinkles
        
        
    elseif strcmpi(G.Opts.Normalization, 'beltrami')
        if G.Opts.Display; fprintf('(beltrami) ... '); end
        
        % beltrami normalization
        D = sparse(1:N, 1:N, DInv, N, N, N);
        K = D*W*D;
        
        DInv = sum(K,2);
        lNonZeroIdxs = find(DInv~=0);
        DInv(lNonZeroIdxs) = DInv(lNonZeroIdxs).^(-1);
        
        V = sparse(1:N, 1:N, DInv, N, N, N);
        T = V*K;
        
    elseif strcmpi(G.Opts.Normalization, 'FokkerPlanck')
        if G.Opts.Display; fprintf('(FokkerPlanck) ... '); end
        
        % beltrami normalization
        D = sparse(1:N, 1:N, sqrt(DInv), N, N, N);
        K = D*W*D;
        
        DInv = sum(K,2);
        lNonZeroIdxs = find(DInv~=0);
        DInv(lNonZeroIdxs) = DInv(lNonZeroIdxs).^(-1);
        
        D = sparse(1:N, 1:N, DInv, N, N, N);
        T = D*K;
        
    elseif strcmpi(G.Opts.Normalization, 'sFokkerPlanck')
        if G.Opts.Display; fprintf('(sFokkerPlanck) ... '); end
        
        % Fokker Planck normalization
        D = sparse(1:N, 1:N, sqrt(DInv), N, N, N);
        K = D*W*D;
        
        DInv = sum(K,2);
        lNonZeroIdxs = find(DInv~=0);
        DInv(lNonZeroIdxs) = DInv(lNonZeroIdxs).^(-1);
        
        D = sparse(1:N, 1:N, sqrt(DInv), N, N, N);
        P = D;
        T = D*K*D;
        
        T = (T+T')/2;           % Iron numerically wrinkles
    else
        fprintf('\nGraphDiffusion:Warning: unknown normalization.');
        return;
    end
    
    if ~strcmpi(G.Opts.Normalization, 'bimarkov')
        if G.Opts.Display; fprintf('%g seconds\n', etime(clock, TIME)); end
    end
    
    if lNumberOfGraphs>1,
        vT{k} = T; vW{k} = W; vP{k} = P;
    end
    
end

if lNumberOfGraphs>1,
    T = vT; W = vW; P = vP;
end

% Save the matrices in the structure G
G.T = T; G.W = W; G.P = P;

if ~G.Opts.DontReturnDistInfo,
    G.DistInfo = DistInfo;
end

% Compute eigenvectors if required
if G.Opts.kEigenVecs>0,
    if G.Opts.Display, fprintf('computing eigenfunctions... '); TIME = clock; end
    [G.EigenVecs, G.EigenVals] = GetEigs(G.T, G.Opts.kEigenVecs, G.P,struct('TakeDiagEigenVals',1));
    if G.Opts.Display, fprintf('%g seconds\n',etime(clock,TIME)); end
end

return;