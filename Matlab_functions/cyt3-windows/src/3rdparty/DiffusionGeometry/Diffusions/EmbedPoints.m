function OutPoints = EmbedPoints(Points, Dimension)
% function OutPoints = EmbedPoints(Points, Dimension)
%
% EMBEDPOINTS trivially embeds a point cloud lying in d-dimension Euclidean
% space into a higher dimensional Euclidean space and then applies a random
% orthogonal matrix to the  set.
%
% In:
%    Points    = DxN matrix representing M points in R^D
%    Dimension = dimension of the embedding space
%
% Out:
%    OutPoints   DimensionxN matrix giving the n embedded points
%
% Dependencies:
%    none
%
% Version History:
%   jcb        1/2006         initial version created
%   jcb        8/2006         points matrix transposed


[d, N] = size(Points);

OutPoints = [Points; zeros(Dimension-d, N)];
OutPoints = RandomOrthogonal(Dimension)'*OutPoints;

function Q = RandomOrthogonal(N)
% RANDOMORTHOGONAL generate uniformly distributed random orthogonal matrixes.
% The distribution is uniformly distributed with respect to Haar measure on
% the group of orthogonal NxN matrices.
A = randn(N);
[Q, R] = qr(A);

