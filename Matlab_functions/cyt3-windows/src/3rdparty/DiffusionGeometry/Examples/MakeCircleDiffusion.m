function T = MakeCircleDiffusion(N)
% function T = MakeCircleDiffusion(N)
%
% MAKECIRCLEDIFFUSION returns the 3-point diffusion operator on the circle.
%
% In:
%    N = number of points for the circle' discretization
%
% Out:
%    T = NxN sparse matrix representing the diffusion
%
% Dependencies:
%    gconv.dll
%
% Version History:
%    jcb       02/06/06       initial version
%

e = 0.5*ones(N,1);
T = spdiags([e e], [-1,1],N,N)

return;
