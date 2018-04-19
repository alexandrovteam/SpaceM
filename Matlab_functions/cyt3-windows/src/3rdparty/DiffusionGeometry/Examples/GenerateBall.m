function X = GenerateBall(n,k,D)

%
% function X = GenerateBall(n,k,D)
% 
% Generates n uniformly distributed points in the k-dimensional unit ball in R^D
%
% IN:
%   n   : number of points
%   k   : dimension of the ball
%   D   : ambient dimension
%
% OUT:
%   X   : n by D matrix of points
%


% Draw points on the unit sphere first (uniform angular)
X       = [randn(n,k),zeros(n,D-k)];
Xnorms  = sqrt(sum(X.^2,2));

for i=1:n
    X(i,:)=X(i,:)/Xnorms(i);
end

% Draw from data distribution for the radial density
beta=betarnd(k,1,[n 1]);

% Push the points in based on the radial distribution
for i=1:n
    X(i,:)=X(i,:)*beta(i);
end

return;
    