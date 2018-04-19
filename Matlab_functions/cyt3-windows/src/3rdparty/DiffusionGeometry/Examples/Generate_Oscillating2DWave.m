function X=Generate_Oscillating2DWave(N, Dim, Opts)

%
% function cX=Generate_Oscillating2DWave(N, Dim, Opts)
%
% Generate_S_Manifold generates a S shape manifold.
% 
% IN:
%    NofPts     : the number of points in the manifold generated
%    [Dim]      : the dimension of the manifold, if Dim=1, a curve, if Dim=2, a surface. default = 2
%    [Opts]     : structure containing the following fields:
%                   [PtsType] : 'mesh': a meshgrid type 
%                               'rand': a random uniform sampling
%                               default = rand
%
% OUT:
%     X: NofPtsxDim array, if PtsType is mesh, not exactly. 
%
%
% SC:
%    MM: 10/2010
%


% Setup parameters
if nargin < 2
    Dim = 2;
end

if nargin < 3
   Opts=[];
end

if Dim > 2
    fprintf('Dim > 2, Dim is modified to 2.\n'); 
    Dim = 2;
end

x = rand(N,1);
x = sort(x);

X = [x,rand(N,1),zeros(N,1)];

X(:,3) = sin( 10*pi*(X(:,1)).^4 );

return;

