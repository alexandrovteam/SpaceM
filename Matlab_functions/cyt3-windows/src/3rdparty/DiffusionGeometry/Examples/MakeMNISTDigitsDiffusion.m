%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Parameters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Which digits to consider
pDigits = [2,3];   
% Number of points
pN = 1000;

%
% Create the data set
fprintf('\n Constructing graph and diffusion map...');
[Data.X,Data.GraphDiffOpts] = GenerateDataSets('BMark_MNIST',...
     struct('NumberOfPoints', pN,'MnistOpts', struct('QueryDigits', pDigits), 'BMarkUnitBall', 1));
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