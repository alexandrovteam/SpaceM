%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Parameters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Intrinsic dimension of the sphere
pK = 2;   
% Number of points
pN = 2000;
% Embedding (ambient) dimension
pD = 3;
% Variance of Gaussian noise to be added
pSigma = 0.01;

%
% Create the data set
fprintf('\n Constructing graph and diffusion map...');
[Data.X,Data.GraphDiffOpts] = GenerateDataSets( 'SwissRoll', ...
    struct('NumberOfPoints',pN,'Dim',pK,'EmbedDim',pD,'NoiseType','Gaussian','NoiseParam',pSigma));
fprintf('done.\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Graph, diffusion construction, eigen computation
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n Constructing graph and diffusion map...');
Data.G = GraphDiffusion(Data.X, 0, Data.GraphDiffOpts);                  
fprintf('done.\n');

Data.Coords = Data.G.EigenVecs;