%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Parameters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Intrinsic dimension of the sphere
pK = 2;
% Number of points
pN = 5000;
% Embedding (ambient) dimension
pD = 3;
% Variance of Gaussian noise to be added
pSigma = 0.01;

%
% Create the data set
fprintf('\n Constructing graph and diffusion map...');
[Data.X,Data.GraphDiffOpts] = GenerateDataSets( 'D-Sphere', ...
    struct('NumberOfPoints',pN,'Dim',pK,'EmbedDim',pD,'NoiseType','Gaussian','NoiseParam',pSigma));
fprintf('done.\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Graph, diffusion construction, eigen computation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n Constructing graph and diffusion map...');
Data.G = GraphDiffusion(Data.X, 0, Data.GraphDiffOpts);                  
fprintf('done.\n');

Data.Coords = Data.G.EigenVecs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Computes Diffusion Wavelets
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fprintf('\n Constructing Diffusion Wavelet Tree...');
% Data.G.WaveletOpts = struct('Wavelets',false,'GS','gsqr_suitesparse');
% Data.Tree = DWPTree (Data.G.T, 20, 1e-6, Data.G.WaveletOpts);
% fprintf('done.\n');


