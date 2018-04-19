function Points = AcceptReject(NPoints, AcceptFcn, Extents)
% function Points = AcceptReject(NPoints, AcceptFcn, [Extents])
%
% ACCEPTREJECT is used to generate a set of uniformly distributed
% points in a arbitrary region R via the "accept-reject" method.  The
% "accept-reject" methods works by choosing uniformly distributed points
% in a bounding rectangle which contains the region R and then accepting
% only those points that happen to fall into R.
%
% ACCEPTREJECT assumes by default that the region falls within the cube
% [-1,1]^d and calls a user provided function with each proposed point.  If
% the point lies in the region R, the function should return 1, and 0
% otherwise.
%
% If the optional Extents = (a_ij) points are chosen within a bounding
% rectangle
%
%    [a_11, a_12] x [a_21, a_22] x ... x [a_d1, a_d2]
%
% instead of over the cube [-1,1]^d.  This can save a significant amount of
% time.
%
% In:
%   NPoints   = number of points to generate
%   AcceptFcn = handle to the reject-accept function
%   Extents   = dx2 array specifying bounds for the bounding rectangle
%
% Out:
%   Points    = set of uniformily distributed points in the desired region
%
% Dependencies:
%   none
%
% Version History:
%   jcb        1/2006         initial version
%
%

T = cputime;

fprintf('AcceptReject.m: generating points ... ');

N = 0;                        % number of points chosen thus far
f = AcceptFcn;                % shortcut for the function
Dimension = nargin(f);        % find out how many arguments the function f has
Pts = 2000;                   % number of test points for each iteration

% allocate space for points
Points     = zeros(Dimension, NPoints);
TempPoints = zeros(Dimension, NPoints);
Accepted = 0;

% this is cheap but it works :)
if exist('Extents')
   a = min(Extents(:,1));
   b = max(Extents(:,2));
else
   a = -1;
   b = 1;
end

fprintf('%05d / %05d', 0, NPoints);

while N < NPoints
   % pick some random points in the cube [a,b]^d
   TempPoints = cell(1, Dimension);
   for j=1:Dimension
      TempPoints{j} = (b-a).*rand(Pts, 1)-b;
   end

   % decide to accept or reject
   Accept = f(TempPoints{:});
   Idxs   = find(Accept);


   K = min(nnz(Accept), NPoints-N);
   Accepted = Accepted+K;

   fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b%05d / %05d', Accepted, NPoints);

   for j=1:Dimension
      Points(j, (N+1):(N+K)) = TempPoints{1,j}(Idxs(1:K), 1);
   end

   N = N+K;


end

fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\bdone (%f ms)\n', (cputime-T)*1000.0);