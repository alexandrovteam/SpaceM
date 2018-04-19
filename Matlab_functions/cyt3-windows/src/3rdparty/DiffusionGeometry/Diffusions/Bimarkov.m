function [T, p] = Bimarkov(K, Options)

%
% function [T, p] = Bimarkov(K, Options)
%
% BIMARKOV computes the bimarkov normalization function p(x) for the
% nonnegative, symemtric  kernel K(x,y) using an iterative scheme.  The
% function p(x) is the unique function s.t.
%
%    diag(1./sqrt(p)*K*diag(1./sqrt(p))
%
% is both row and column stochastic.  Note that a bimarkov kernel
% necessarily has L^2 norm 1.
%
% In:
%    K       = an NxN matrix speciying a nonnegative, symmetric kernel with
%              nonzero row sums
%    Options = a structure containing zero or more of the following fields
%
%      MaxIters : maximum number of iterations (default = 100)
%      AbsError : BIMARKOV stops when each row sum is within AbsoluteError
%                      of 1.00 (default = 1e-5)
%      Quiet    : when true, all output other than error messages will be
%                 suppressed
%
% Out:
%    T       = bimarkov normalized kernel
%    p       = column vector giving the bimarkov normalization function
%
% Dependencies:
%    none
%
% Version History:
%   jcb        12/2006        initial version
%   mm         4/2007         some small bug fixes
%
% (c) Copyright Yale University, 2006, James C Bremer Jr.
% (c) Copyright Duke University, 2007, Mauro Maggioni
%


TIME = cputime;

% initialize output
T = [];
p = [];

% nothing in, nothing out
if isempty(K)
   return;
end

% process input
if size(K,1)~=size(K,2)
   fprintf('Bimarkov.m: kernel must be NxN\n');
   return;
end

if nargin<2,
   Options = [];
end

if ~isfield(Options, 'MaxIters')
   Options.MaxIters = 100;
end

if ~isfield(Options, 'AbsError')
   Options.AbsError = 1e-5;
end
if ~isfield(Options, 'Quiet')
   Options.Quiet = 0;
end


N = size(K,1);
MaxIters = Options.MaxIters;
AbsError = Options.AbsError;


if ~Options.Quiet
   fprintf('Bimarkov: normalizing %dx%d kernel ... ', N,N);
end

% initialize
p = ones(N, 1);

E = cell(MaxIters, 1);

% iterate
for iter=1:MaxIters

   S = full(sum(K,2));
   err = max(abs(1.0-max(S)), abs(1.0 - min(S)));

   if ~Options.Quiet
      fprintf('iter %03d  error: 1e%1.4f', iter, log10(err));
   end

   if err < AbsError
      if ~Options.Quiet
         if log10(err) < 0
            fprintf('\b');
         end
         fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b');
      end
      break;
   end

   D = sparse(1:N, 1:N, 1./sqrt(S), N, N, N);
   p = S.*p;
   K = D*K*D;

   if ~Options.Quiet && log10(err) < 0
      fprintf('\b');
   end

   if ~Options.Quiet
      fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b');
   end
end

if ~Options.Quiet
   fprintf('\b done (%4.2f ms, iters = %d, err = %g) \n', (cputime-TIME)*1000.0, iter, err);
end

% iron out numerical errors
T = (K'+K)/2;


return;