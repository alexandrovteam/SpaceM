function [X,GraphDiffOpts,NetsOpts,Labels,Coords,GWTOpts] = GenerateDataSets( cDataName, cOpts )

%
% function [X,GraphDiffOpts,NetsOpts,Labels,Coords] = GenerateDataSets( cDataName, cOpts )
%
% IN:
%   cDataName   : name of the data set to be generated. Known names:
%                   'D-Cube'        : D-dimesional cube
%                   'D-Ball'        : D-dimensional unit ball. Quite inefficient in high-dimensions.,
%                   'D-Sphere'      : D-dimensional sphere in D dimensions
%                   'D-FlatTorus'   : D-dimensional flat torus in 2D dimensions
%                   'D-Gaussian'    : D-dimensional Gaussian distribution One may impose a standard deviation via [GaussianStdDev] in cOpts.
%                                       Default standard deviation = 1
%                   'Paraboloid'    : 1/2 x^T H x, for a positive semidefinite Hessian matrix H
%                   'S-Manifold'    : S-shape manifold (1d-curve, 2d-surface)
%                   'Oscillating2DWave' : a plane that oscillates faster and faster
%                   'SwissRoll'     : Swiss-roll-shape manifold (1d-curve, 2d-surface)
%                   'ImagePatches'  : patches from an image, from a file name cOpts.ImageFileName, and patch size cOpts.PatchSize
%                   'FromImage'     : data from a two-dimensional image, from a file name cOpts.ImageFileName.
%                   'MeyerStaircase': Meyer's staircase example
%                   'MeyerStaircase-d': Meyer's staircase example, generalized to d dimensions
%                   'Planes'        : planes of various dimensions in various positions (see cOpts.PlanesOpts)
%                   'SpiralAndPlane': spiral and plane data set from the paper by Haro et al. on Intrinsic Dimensionality Estimation with Translated Poisson Processes
%                   'SphereAndLine' : a 2d sphere intersecting a line
%                   'TwoLinesAndAPlane' : 200 points on each of two lines, and 400 points on a portion of a plane, interesecting the two intersecting lines.R
%
%                 The following names of data wil not be generated, instead will be loaded. It will be assumed that those data set are in
%                 the same directory or in the search path. More options are possible, read the headers of Generate_S_Manifold.m,
%                 Generate_Swissroll.m, Generate_MNIST.m. In those case, Dim=EmbedDim= 241 for the Chapelle's data, 28^2 for the MNIST
%
%                   'BMark_g241c'   : g241c in the Chapelle's benchmark set
%                   'BMark_g241d'   : g241d    "
%                   'BMark_Digit1'  : Digit1   "
%                   'BMark_USPS'    : USPS     "
%                   'BMark_COIL'    : COIL     "
%                   'BMark_BCI'     : BCI      "
%                   'BMark_Text'    : Text     "
%                   'BMark_MNIST'   : the MNIST database of handwritten digits
%                   'CBCLFaces1'   : CBCL Faces DataBase
%                   'CBCLFaces2'   : another set from the CBCL Faces DataBase
%                   'HandVideo'    :
%                   'FaceVideo'    :
%                   'Isomapfaces'  :
%                   'ScienceNews'  : 1161 Science News articles, prepared by J. Solka. It normalizes to unit L^2 norm (correlation distance).
%                   'Cosine'       : cosines of frequency 2\pi k, k ranging from 1 to cOpts.NumberOfPoints, sampled at cOpts.Dim equispaced points in [0,1]
%                   'Signal1D_1'
%                   'NaturalImagePatches' : a data set of about 75001 patches of size 16x16 pixels from multiple natural images
%                   'SquareTranslations' : translates the contents of a data cube (with wrap-arounds).
%                   'STL10train'      : STL-10 train data set.
%                   'STL10test'       : STL-10 test data set.
%                   'STL10unlabeled'  : STL-10 unlabeled data set.
%                   'MM_Faces'        : pictures of faces
%
%   cOpts       : structure of options:
%                   [NumberOfPoints]  : number of data points. Default: 1000.
%                   [Dim]             : dimension of the data.
%                   [EmbedDim]        : embedding dimension. Default: equal to data set dimension.
%                   [NoiseType]       : 'uniform' or 'Gaussian'. Default: Gaussian.
%                   [NoiseLocation]   : 'tangent', 'normal', or 'entire'. Default: tangent
%                                       'tangent': noise on the data set
%                                       'normal' : noise on the normal space to the data set in the embedded space
%                                       'entire' : noise on the entire embedded space
%                   [NoiseParam]      : size of the noise (dilation factor). Default: 0.
%                   [GaussianMean]    : mean for Gaussian dist, only for 'D-Gaussian' case. Default: 0.
%                   [GaussianStdDev]  : std deviation for Gaussian dist, only for 'D-Gaussian' case. Default: 1
%                   [BMarkUnitBall]   : Normalize data onto the unit ball by dividing by its maximum norm. Only for real world data sets. Default: 0.
%                   [MnistOpts]       : Options for the MNIST case. See the header of Generate_MNIST.m
%                   [ImageFileName]   : if cDataName is 'ImagePatches' or 'FromImage', it loads this file. Black dots in the image become 2-D points.
%                   [PatchSize]       : if cDataName is 'ImagePatches', extracts patches of this size from the image.
%                   [DownsampleImage] : if cDataName is 'ImagePatches', this is the downsampling factor for the image before extracting the patches.
%                   [AutotuneScales]  : compute set of scales for the data set to go in NetsOpts. Default: true if NetsOpts is an output.
%                   [MeyerStepWidth]  : how wide a step in the MeyerStairCase is. Default: 10.
%                   [MeyerStepNorm]   : p for the p-norm to be used in the generation of teh MeyerStairCase. Default: 2.
%                   [DilationFactors] : cOpts.EmebedDim vector of dilation factors, one per direction, applied after noise is added.
%                   [RandomDilationFactors] : a vector of dilation factors, from which uniform draws are taken to dilate each axis. This is done after DilationFactors are applied. Default: [1].
%                   [RandomScaling]   : a vector of global isotropic scaling factors for all the axes, from which a uniform draw is taken. Defult: [1].
%                   [SVDProject]      : Project data, after noise added, onto this number of SVD vectors. Default: 0 (no projection).
%                   [PlanesOpts]      : structure of options for planes:
%                           Number      : number of planes
%                           Dim         : dimension of planes (either a vector of length Number or a singleton)
%                           Origins     : Number*EmbedDim matrix of centers of the planes
%                           Dilations   : how much to dilate each plane
%                           NumberOfPoints : vector of length Number indicating the number of points on the planes. Should sum to cOpts.NumberOfPoints.
%                   [ParaboloidOpts]:
%                           H           : Hessian matrix
%                   [RandomProjectionDim]: If >0, applies a random projection in RandomProjectionDim dimensions. This is done after everything else (SVDProject, addition of noise, etc...). Default: 0.
%
% OUT:
%   X             : Dim by NumberOfPoints matrix of data points. It is set as global variable.
%   GraphDiffOpts : structure of options for GraphDiffusion.
%   NetsOpts      : structure of options for FastRoughMultiscaleSVD.
%   Labels        : NumberOfPoints column vector of labels for the data points. It is set as global variable.
%   Coords        : Coordinates for X. For example if SVDProject is specified, Coords is the basis (each basis element in a column) onto which X was projected.
%
%
%
% USES:
%   DataPreprocess.m (BMarkUnitBall option for all benckmarksets)
%   Generate_Gaussian.m (D-Gaussian only),
%   Generate_S_Manifold.m (S-manifold only),
%   Generate_Swissroll.m (Swiss roll case only)
%   Generate_MNIST.m (MNIST handwritten digits only)
%   Data: SSL,set=*,data.mat (Chapelle's benchmark sets only)
%         TrainImages.mat, TrainImageLabels.mat (MNIST data only)
%
% Copyright (c)
% Duke University, 2008
% Mauro Maggioni
%

Coords = [];
Labels = [];
GraphDiffOpts = [];
NetsOpts = [];

AutoScalesMinNPerBin = 10;

if nargin<2,
    cOpts = [];
end
if ~isfield(cOpts,'NumberOfPoints'),        cOpts.NumberOfPoints = 1000;        end
if ~isfield(cOpts,'Dim'),                   cOpts.Dim = 0;                      end
if ~isfield(cOpts,'EmbedDim'),              cOpts.EmbedDim = cOpts.Dim;         end
if ~isfield(cOpts,'NoiseType'),             cOpts.NoiseType = 'Gaussian';       end
if ~isfield(cOpts,'NoiseParam'),            cOpts.NoiseParam = 0;               end
if ~isfield(cOpts,'MeyerStepWidth'),        cOpts.MeyerStepWidth = 10;          end
if ~isfield(cOpts,'MeyerStepNorm'),         cOpts.MeyerStepNorm = 2;            end
if ~isfield(cOpts,'DilationFactors'),       cOpts.DilationFactors = [];         end
if ~isfield(cOpts,'SVDProject'),            cOpts.SVDProject = 0;               end
if ~isfield(cOpts,'RandomProjectionDim'),   cOpts.RandomProjectionDim = 0;      end
if ~isfield(cOpts,'RandomDilationFactors'), cOpts.RandomDilationFactors = [];   end
if ~isfield(cOpts,'RandomScaling'),         cOpts.RandomScaling = [];           end
if ~isfield(cOpts,'AutotuneScales'),        cOpts.AutotuneScales = true;        end

% Set default options for constructing the graph
GraphDiffOpts = struct( ...
    'Normalization','smarkov', ...
    'Epsilon',1, ...
    'kNN', 100, ...
    'kNNAutotune', 20, ...
    'kEigenVecs', 35, ...
    'Symmetrization', 'W+Wt', ...
    'DontReturnDistInfo', 1 );
% Set default options for constructing the multiscale svd nets
NetsOpts = struct( ...
    'Epsilon',1, ...
    'Delta',[0.05:0.025:1], ...
    'NumberOfScales',length([0.05:0.025:1]), ...
    'SVDDim',cOpts.EmbedDim, ...
    'ApproxSVD',0, ...
    'NormalizedSVD',1, ...
    'DecreasingNets',0, ...
    'DownSample',1, ...
    'KeepV',0);
% Set default options for constructing geometric wavelets
GWTOpts = struct ( ...
    'AmbientDimension', cOpts.EmbedDim, ...
    'ManifoldDimension', cOpts.Dim, ...
    'errorType', 'relative', ...
    'precision', 1e-2, ...
    'threshold0', 0.5, ...
    'threshold1', 1e-4, ...
    'threshold2', 1e-4, ...
    'sparsifying', false, ...
    'splitting', false, ...
    'pruning', 1 ...
    );

switch lower(cDataName)
    case {'d-cube'}
        Xo = rand(cOpts.Dim,cOpts.NumberOfPoints);
    case {'d-ball'}
        Xo = GenerateBall(cOpts.NumberOfPoints,cOpts.Dim,cOpts.Dim);
        Xo = Xo';
    case {'d-sphere'}
        Xo = randn(cOpts.NumberOfPoints,cOpts.Dim+1);
        for k = 1:cOpts.NumberOfPoints,
            Xo(k,:) = Xo(k,:)/norm(Xo(k,:));
        end
        Xo = Xo';
    case {'sphereandline'}
        lOpts = cOpts;
        lOpts.NumberOfPoints = round(2/3*cOpts.NumberOfPoints);
        lOpts.Dim = 2;
        lOpts.EmbedDim = 3;
        [X_1,GraphDiffOpts,NetsOpts] = GenerateDataSets( 'D-Sphere',lOpts );
        lOpts.NumberOfPoints = round(1/3*cOpts.NumberOfPoints);
        lOpts.Dim = 1;
        lOpts.EmbedDim = 3;
        X_2 = GenerateDataSets( 'D-Cube', lOpts );
        Xo = [X_1,X_2];
        Labels = [ones(size(X_1,1),1);2*ones(size(X_2,1),1)];
    case {'twolinesandaplane'}
        load TwoLinesAndAPlane;
        Xo = X'; clear X;
    case {'d-flattorus'}
        Xo   = zeros(cOpts.NumberOfPoints,2*cOpts.Dim);
        lTmp = rand(cOpts.NumberOfPoints,cOpts.Dim);
        for d = 1:cOpts.Dim,
            Xo(:,2*d-1) = cos(2*pi*lTmp(:,d));
            Xo(:,2*d)   = sin(2*pi*lTmp(:,d));
        end
        Xo = Xo';
    case {'d-gaussian'}
        if isfield(cOpts,'GaussianStdDev'),
            if length(cOpts.GaussianStdDev)>1
                Xo = Generate_Gaussian(cOpts.NumberOfPoints, cOpts.Dim, cOpts.GaussianStdDev);
            else
                Xo = cOpts.GaussianStdDev*randn(cOpts.NumberOfPoints, cOpts.Dim);
            end
        else
            Xo = randn(cOpts.NumberOfPoints, cOpts.Dim);
            cOpts.GaussianStdDev = 1;
        end
        
        if isfield(cOpts,'GaussianMean'),
            Xo = bsxfun(@plus,Xo,cOpts.GaussianMean);
        end
        Xo = Xo';
    case {'paraboloid'}
        Xo = [2*rand(cOpts.NumberOfPoints,cOpts.Dim)-1,zeros(cOpts.NumberOfPoints,1)];
        if ~isfield(cOpts,'ParaboloidOpts'),
            cOpts.ParaboloidOpts.H = eye(cOpts.Dim);
        end
        for k = 1:cOpts.NumberOfPoints,
            Xo(k,cOpts.Dim+1) = 1/2*Xo(k,1:cOpts.Dim)*cOpts.ParaboloidOpts.H*Xo(k,1:cOpts.Dim)';
        end
        Xo = Xo';
    case {'s-manifold'}
        Xo = Generate_S_Manifold(cOpts.NumberOfPoints,cOpts.Dim);
        Xo = Xo';
    case {'oscillating2dwave'}
        Xo = Generate_Oscillating2DWave(cOpts.NumberOfPoints,cOpts.Dim);
        Xo = Xo';
    case {'swissroll'}
        Xo = Generate_Swissroll(cOpts.NumberOfPoints,cOpts.Dim);
        Xo = Xo';        
    case {'spiralandplane'}
        Xo = load('Xn3_EspPlano_9oct06.mat');
        Xo = Xo.X3n;
    case {'meyerstaircase-d'}        
        nptsperdim = round(cOpts.EmbedDim^(1/cOpts.Dim));
        ptsperdim  = linspace(-1,1,nptsperdim);
        
        if ~isfield(cOpts,'width'),
            cOpts.width = 1/(nptsperdim/3);
        end
        if ~isfield(cOpts,'MeyerStepNorm')
            cOpts.MeyerStepNorm = 2;
        end
        % Create grid in [-1,1]^D with cOpts.NumberOfPoints points
        ntauperdim = round(cOpts.NumberOfPoints^(1/cOpts.Dim));
        tauperdim  = linspace(-1,1,ntauperdim);        
        for k = 1:cOpts.Dim,
            str_eval = '[';
            str_eval = [str_eval,sprintf('Xs_%d',k),']=rand('];
            for j = 1:cOpts.Dim,
                str_eval = [str_eval,num2str(ntauperdim),','];
            end
            str_eval(end) = ',';
            str_eval = [str_eval,'1)*2-1;'];
            eval(str_eval);
        end
                
        % Construct Gaussians
        Xo          = zeros(prod(size(Xs_1)),nptsperdim^cOpts.Dim);
        Coords      = zeros(prod(size(Xs_1)),cOpts.Dim);
        idx = 1;
        for i = 1:prod(size(Xs_1)),
            str_eval   = 'Coords(idx,:) = [';
            for k = 1:cOpts.Dim,
                str_eval   = [str_eval sprintf('Xs_%d(i),',k)];
            end
            str_eval(end)  = [];
            str_eval       = [str_eval '];'];
            eval(str_eval);
            tmp=CreateGriddedGaussian(cOpts.Dim,ptsperdim,Coords(idx,:),cOpts.width, cOpts.MeyerStepNorm);
            Xo(idx,:)       = tmp(:)/max(abs(tmp(:)));
            idx             = idx + 1;
        end
        Labels = Coords;
        AutoScalesMinNPerBin = 10;     
        Xo = Xo';        
    case {'meyerstaircase'}
        Xo = zeros(cOpts.NumberOfPoints,cOpts.Dim);
        for k = 1:cOpts.NumberOfPoints,
            Xo(k,k:min([k+cOpts.MeyerStepWidth,cOpts.Dim])) = 1;
            Xo(k,:) = Xo(k,:)/norm(Xo(k,:));
            Xnew(k,:) = conv(Xo(k,:),conv(conv(Xo(k,:),Xo(k,:),'full'),Xo(k,:),'full'),'full');
        end
        Xo = Xnew;
        Labels = (1:cOpts.NumberOfPoints)';
        AutoScalesMinNPerBin = 10;
        GraphDiffOpts.kNN = 5;
        GraphDiffOpts.kNNAutotune = 3;
        Xo = Xo';        
    case {'bmark_g241c'}
        load ('SSL,set=5,data.mat');
        cOpts.Dim=size(X, 2); cOpts.EmbedDim=size(X, 2);
        if isfield(cOpts,'BMarkUnitBall') && cOpts.BMarkUnitBall,
            X=DataPreprocess(X, struct('UnitBall', 1));
        end
        Xo=X'; clear X y;
    case {'bmark_g241d'}
        load ('SSL,set=7,data.mat');
        cOpts.Dim=size(X, 2); cOpts.EmbedDim=size(X, 2);
        if isfield(cOpts,'BMarkUnitBall') && cOpts.BMarkUnitBall,
            X=DataPreprocess(X, struct('UnitBall', 1));
        end
        Xo=X'; clear X y;
    case {'bmark_digit1'}
        load ('SSL,set=1,data.mat');
        cOpts.Dim=size(X, 2); cOpts.EmbedDim=size(X, 2);
        if isfield(cOpts,'BMarkUnitBall') && cOpts.BMarkUnitBall,
            X=DataPreprocess(X, struct('UnitBall', 1));
        end
        Xo=X'; clear X y;
    case {'bmark_usps'}
        load ('SSL,set=2,data.mat');
        cOpts.Dim=size(X, 2); cOpts.EmbedDim=size(X, 2);
        if isfield(cOpts,'BMarkUnitBall') && cOpts.BMarkUnitBall,
            X=DataPreprocess(X, struct('UnitBall', 1));
        end
        Xo=X'; clear X y;
    case {'bmark_coil'}
        load ('SSL,set=6,data.mat');
        cOpts.Dim=size(X, 2); cOpts.EmbedDim=size(X, 2);
        if isfield(cOpts,'BMarkUnitBall') && cOpts.BMarkUnitBall,
            X=DataPreprocess(X, struct('UnitBall', 1));
        end
        Xo=X'; clear X y;
    case {'bmark_bci'}
        load ('SSL,set=4,data.mat');
        cOpts.Dim=size(X, 2); cOpts.EmbedDim=size(X, 2);
        if isfield(cOpts,'BMarkUnitBall') && cOpts.BMarkUnitBall,
            X=DataPreprocess(X, struct('UnitBall', 1));
        end
        Xo=X'; clear X y;
    case {'bmark_text'}
        load ('SSL,set=9,data.mat');
        X=full(X);
        cOpts.Dim=size(X, 2); cOpts.EmbedDim=size(X, 2);
        if isfield(cOpts,'BMarkUnitBall') && cOpts.BMarkUnitBall,
            X=DataPreprocess(X, struct('UnitBall', 1));
        end
        Xo=X'; clear X y;
    case {'bmark_mnist'}
        if ~isfield(cOpts, 'MnistOpts')
            [Xo,Labels]=Generate_MNIST(cOpts.NumberOfPoints);
        else
            [Xo,Labels]=Generate_MNIST(cOpts.NumberOfPoints, cOpts.MnistOpts);
        end
        if isfield(cOpts,'BMarkUnitBall') && cOpts.BMarkUnitBall,
            Xo=DataPreprocess(Xo, struct('UnitBall', 1));
        end
        cOpts.Dim=size(Xo, 2); cOpts.EmbedDim=size(Xo, 2);
        Xo = Xo';
    case {'fromimage'}
        try
            % Load the file containing the description of the rooms environment
            lImage = imread(cOpts.ImageFileName);
        catch
            fprintf('\n Error: could not load the description file %s!',cOpts.ImageFileName);
            return;
        end
        
        lSumImage = sum(lImage,3);
        % Find the points in the image
        [Xo(1,:),Xo(2,:)] = find(lSumImage==0);
        
        cOpts.Dim=size(Xo, 2); cOpts.EmbedDim=size(Xo, 2);
        Xo = Xo';
        
        lRoughDiam = max(max(Xo))-min(min(Xo));
        lDeltas = 2.^linspace(log2(lRoughDiam/(2^10)),log2(lRoughDiam/2),15);
        
        GraphDiffOpts = struct( ...
            'Normalization','smarkov', ...
            'Epsilon',1, ...
            'kNN', 20, ...
            'kNNAutotune', 10, ...
            'kEigenVecs', 35, ...
            'Symmetrization', 'W+Wt', ...
            'DontReturnDistInfo', 1 );
        % Set default options for constructing the multiscale svd nets
        NetsOpts = struct( ...
            'Epsilon',1, ...
            'Delta',lDeltas, ...
            'NumberOfScales',length(lDeltas), ...
            'SVDDim',cOpts.EmbedDim, ...
            'ApproxSVD',0, ...
            'NormalizedSVD',1, ...
            'DecreasingNets',0, ...
            'DownSample',1, ...
            'KeepV',0);        
    case {'imagepatches'}
        try
            % Load the file containing the description of the rooms environment
            lImage = imread(cOpts.ImageFileName);
        catch
            fprintf('\n Error: could not load the description file %s!',cOpts.ImageFileName);
            return;
        end
        
        lImage=double(lImage);
        lImage=lImage/max(max(lImage));
        
        %% Downsample image
        if isfield(cOpts,'DownsampleImage') && ~isempty(cOpts.DownsampleImage),
            lImage = lImage(1:cOpts.DownsampleImage:size(lImage,1),1:cOpts.DownsampleImage:size(lImage,2));
        end
        
        %% Construct data set
        if ~isfield(cOpts,'PatchSize') || isempty(cOpts.PatchSize),
            cOpts.PatchSize = 16;
        end
        Xo = FilterGraph(lImage,'ptch',cOpts.PatchSize);
        Xo = Xo';
        
        cOpts.Dim=size(Xo, 2); cOpts.EmbedDim=size(Xo, 2);
        
        GraphDiffOpts = struct( ...
            'Normalization','smarkov', ...
            'Epsilon',1, ...
            'kNN', 20, ...
            'kNNAutotune', 10, ...
            'kEigenVecs', 35, ...
            'Symmetrization', 'W+Wt', ...
            'DontReturnDistInfo', 1 );
        % Set default options for constructing the multiscale svd nets
        NetsOpts = struct( ...
            'Epsilon',1, ...
            'SVDDim',cOpts.EmbedDim, ...
            'ApproxSVD',0, ...
            'NormalizedSVD',1, ...
            'DecreasingNets',0, ...
            'DownSample',1, ...
            'KeepV',0);
    case {'planes'}
        if length(cOpts.PlanesOpts.Dim)==1,
            cOpts.PlanesOpts.Dim = cOpts.PlanesOpts.Dim*ones(cOpts.PlanesOpts.Number);
        end
        if ~isfield(cOpts.PlanesOpts,'Dilations'),
            cOpts.PlanesOpts.Dilations = 1;
        end
        if length(cOpts.PlanesOpts.Dilations)==1,
            cOpts.PlanesOpts.Dilations = cOpts.PlanesOpts.Dilations*ones(cOpts.PlanesOpts.Number);
        end
        
        Xo      = zeros(cOpts.NumberOfPoints,cOpts.EmbedDim);
        Labels  = zeros(cOpts.NumberOfPoints,1);
        idx     = 1;
        
        for k = 1:cOpts.PlanesOpts.Number
            Xtmp = randn(cOpts.EmbedDim,cOpts.PlanesOpts.Dim(k));
            [Qo,Ro] = qr(Xtmp,0);
            Xo(idx:idx+cOpts.PlanesOpts.NumberOfPoints(k)-1,:) = (cOpts.PlanesOpts.Dilations(k)*(rand(cOpts.PlanesOpts.NumberOfPoints(k),cOpts.PlanesOpts.Dim(k))-0.5))*Qo'+repmat(cOpts.PlanesOpts.Origins(k,:),cOpts.PlanesOpts.NumberOfPoints(k),1);
            Labels(idx:idx+cOpts.PlanesOpts.NumberOfPoints(k)-1,:) = k;
            idx = idx+cOpts.PlanesOpts.NumberOfPoints(k);
        end
        Xo = Xo';

    case {'facevideo'}
        load ('FaceVideo.mat');
        cOpts.EmbedDim=size(X, 2); cOpts.Dim=size(X, 2);
        if (isfield(cOpts,'RestrictTrainingSet')),
            X=X(1:cOpts.RestrictTrainingSet,:);
        end
        if isfield(cOpts,'BMarkUnitBall') && cOpts.BMarkUnitBall,
            X=DataPreprocess(X, struct('UnitBall', 1));
        end
        Xo=X'; clear X;
    case {'isomapfaces'}
        load ('IsomapFaces.mat');
        cOpts.EmbedDim=size(X, 2);         cOpts.Dim=size(X, 2);
        if (isfield(cOpts,'RestrictTrainingSet')),
            X=X(1:cOpts.RestrictTrainingSet,:);
        end
        if isfield(cOpts,'BMarkUnitBall') && cOpts.BMarkUnitBall,
            X=DataPreprocess(X, struct('UnitBall', 1));
        end
        Xo=X'; clear X;
    case {'handvideo'}
        load ('HandVideo.mat');
        cOpts.EmbedDim=size(X, 2);        cOpts.Dim=size(X, 2);
        if (isfield(cOpts,'RestrictTrainingSet')),
            X=X(1:cOpts.RestrictTrainingSet,:);
        end
        if isfield(cOpts,'BMarkUnitBall') && cOpts.BMarkUnitBall,
            X=DataPreprocess(X, struct('UnitBall', 1));
        end
        Xo=X'; clear X;
    case {'cbclfaces1'}
        load ('CBCLFaces1.mat');
        cOpts.EmbedDim=size(X, 2);        cOpts.Dim=size(X, 2);
        if (isfield(cOpts,'RestrictTrainingSet')),
            X=X(1:cOpts.RestrictTrainingSet,:);
            X=X(1:cOpts.NumberOfPoints,:);
        end
        if isfield(cOpts,'BMarkUnitBall') && cOpts.BMarkUnitBall,
            X=DataPreprocess(X, struct('UnitBall', 1));
        end
        Xo=X'; clear X;
    case {'cbclfaces2'}
        load ('CBCLFaces2.mat');
        cOpts.EmbedDim=size(X, 2);        cOpts.Dim=size(X, 2);
        if (isfield(cOpts,'RestrictTrainingSet')),
            X=X(1:cOpts.RestrictTrainingSet,:);
        end
        if isfield(cOpts,'BMarkUnitBall') && cOpts.BMarkUnitBall,
            X=DataPreprocess(X, struct('UnitBall', 1));
        end
        Xo=X'; clear X;
    case {'sciencenews'}
        load ('ScienceNews.mat');
        cOpts.EmbedDim=size(X, 2);
        cOpts.Dim=size(X, 2);
        if (isfield(cOpts,'RestrictTrainingSet')),
            X=X(1:cOpts.RestrictTrainingSet,:);
        end
        X=DataPreprocess(X, struct('UnitSphere', 1));
        Xo=X'; clear X;
        Labels = classes(:,1);
    case {'cosine'}
        Xo = zeros(cOpts.NumberOfPoints,cOpts.Dim);
        xvalues = 2*pi*linspace(0,1,cOpts.Dim);
        for k = 1:cOpts.NumberOfPoints,
            Xo(k,:) = cos(k*xvalues);
            Xo(k,:) = Xo(k,:)/norm(Xo(k,:));
        end
        Xo = Xo';
        
    case {'signal1d_1'}
        Xo = zeros(cOpts.NumberOfPoints,cOpts.Dim);
        xvalues = 2*pi*linspace(0,1,cOpts.NumberOfPoints);
        yvalues(find(xvalues<pi)) = sin(2*pi*xvalues(find(xvalues<pi)));
        yvalues(find(xvalues>=pi)) = 1+sin(2*pi*xvalues(find(xvalues>=pi)));
        Xo = PatchifySignal(yvalues,cOpts.Dim,1);
        Xo = Xo';        
    case {'naturalimagepatches'}
        load('NaturalImagePatches');
        Xo = Y; clear Y;
    case {'squaretranslations'}
        lCubeSide = round((cOpts.NumberOfPoints)^(1/2));
        Xo = zeros(lCubeSide^2,lCubeSide,lCubeSide);
        idx = 1;
        for i = 1:lCubeSide,
            for j = 1:lCubeSide,
                Xo(idx,i:min([i+10,lCubeSide]),j:min([j+10,lCubeSide]))=1;
                idx = idx+1;
            end
        end
        Xo = reshape(Xo,[lCubeSide^2,lCubeSide^2]);
        Xo = Xo';        
    case {'stl10train'}
        load('STL10-train');
        Xo = X;
        clear X;
        Xo = double(Xo);
        Labels = y;
    case {'stl10test'}
        load('STL10-test');
        Xo = X;
        clear X;
        Xo = double(Xo);
        Labels = y;
    case {'stl10unlabeled'}
        load('STL10-unlabeled');
        Xo = X;
        clear X;
        Xo = double(Xo);
        Labels = y;
    case {'stl10traing'}
        load('STL10-traing');
        Xo = X;
        clear X;
        Xo = double(Xo);
        Labels = y;
    case {'stl10testg'}
        load('STL10-testg');
        Xo = X;
        clear X;
        Xo = double(Xo);
        Labels = y;
    case {'stl10unlabeledg'}
        load('STL10-unlabeledg');
        Xo = X;
        clear X;
        Xo = double(Xo);
        Labels = y;
    case {'mm_faces'}
        img_size = floor(sqrt(cOpts.Dim));
        [Xo,Labels,LabelNames]  = CreateLabelledRandomFaces( img_size, cOpts.NumberOfPoints ); 
        Xo=double(Xo);
    otherwise
        fprintf('GenerateDataSets: unknown DataName (%s).\n',cDataName)
        return;
end

% Project onto SVD if requested
if cOpts.SVDProject>0,
    meanXo = mean(Xo,2);
    bsxfun(@minus,Xo,meanXo);
    [lU,~,~] = svd(Xo,0);
    Coords = lU(:,1:cOpts.SVDProject);
    Xo = Coords'*Xo;
    NetsOpts.Delta = NetsOpts.Delta*sqrt(cOpts.SVDProject/cOpts.Dim);
    if cOpts.EmbedDim == cOpts.Dim,
        cOpts.EmbedDim = cOpts.SVDProject;
    end
end

cOpts.SVDDim = min([cOpts.Dim*2,cOpts.EmbedDim]);

% Apply deterministic DilationFactors
if ~isempty(cOpts.DilationFactors),
    for k = 1:size(X,2),
        Xo(k,:) = Xo(k,:)*cOpts.DilationFactors(k);
    end
end

% Apply RandomDilationFactors
if ~isempty(cOpts.RandomDilationFactors),
    lRandDilations = cOpts.RandomDilationFactors(randi(length(cOpts.RandomDilationFactors),size(Xo,2),1));
    for k = 1:size(Xo,2),
        Xo(k,:) = Xo(k,:)*lRandDilations(k);
    end
    NetsOpts.RandomDilationFactorsRealization = lRandDilations;
end

% Embed in high dimensions
if isnumeric(cOpts.EmbedDim)
    if cOpts.EmbedDim-size(Xo,1)>0
        Xo = [Xo;zeros(cOpts.EmbedDim-size(Xo,1),size(Xo,2))];
    end
end

% Add noise
if cOpts.NoiseParam > 0
    switch lower(cOpts.NoiseType)
        case {'uniform'}
            X = Xo + cOpts.NoiseParam*(rand(size(Xo))-0.5);
        case {'gaussian'}
            X = Xo + cOpts.NoiseParam*randn(size(Xo));
        otherwise
            fprintf('GenerateDataSets: unknown NoiseType (%s).\n',cOpts.NoiseType);
    end
else
    X=Xo;
end

clear Xo;

% Random projection
if cOpts.RandomProjectionDim >0,
    [Qrnd,~] = qr(randn(cOpts.RandomProjectionDim,cOpts.EmbedDim));
    X = Qrnd*X;
end

% Random scaling
if ~isempty(cOpts.RandomScaling),
    lRandDilation = cOpts.RandomScaling(randi(length(cOpts.RandomScaling),1,1));
    X = lRandDilation*X;
    NetsOpts.RandomScalingRealization = lRandDilation;
end

if nargout>2
    % Autotune the scales
    if cOpts.AutotuneScales,
        NetsOpts.Delta = AutoScales(X,struct('MinNperBin',min([min(size(X,1)),AutoScalesMinNPerBin]),'StatType','Min'));
        NetsOpts.NumberOfScales = length(NetsOpts.Delta);
    end
end

return;
