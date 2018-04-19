function p = cartprod(a,b)
% function p = cartprod(a,b)
%
% CARTPROD forms the cartesian product of the vectors a and b.  That is, the
% length(a)*length(b) x 2 matrix p such that
%
%    p(i,j) = [a(ceil(i/length(a)) b(i%length(a)].
%
% In:
%   a = vector of length k
%   b = vector of length l
%
% Out:
%   p = (k*l)x2 cartesian product matrix
%
% Version History:
%

if size(a, 2) > 1
    a=a';
end

if size(b,2) > 1
    b=b';
end

% initialize the output matrix
p = zeros(length(a)*length(b), 2);

% fill in the matrix by varying the values for one set
index = 1;
blocklen = length(b);

for j=1:length(a)
    p(index:(index+blocklen-1), 1) = a(j)*ones(blocklen, 1);
    p(index:(index+blocklen-1), 2) = b;
    index = index + blocklen;
end

