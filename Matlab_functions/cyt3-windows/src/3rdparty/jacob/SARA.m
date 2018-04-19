% Statistical Analysis of Response Amplitude
% -----------------------------------------------------------------------
% Score response distribution y from unperturbed distribution x,
%  using niter randomized values in null distribution
% 
% [score,pval,obs,mdiff,null] = SARA( x, y, niter )

function [score,pval,obs,mdiff,null] = SARA( x, y, niter )
if nargin < 3
	niter = 1000;
end

% initialize null
null = nan(niter,1);

% pool data
xy = [x; y];
n = length(xy);
nn = length(x);

% use fixed equipartition of space
xedges = linspace(min(xy),max(xy),100);

% observed statistic
obs = emd_calc( x, y, xedges );

% build null
stop = false;
M = 0;
ii = 1;
while ~stop && ii < niter+1
    rix = randperm(n);
    null(ii,1) = emd_calc( xy(rix(1:nn)), xy(rix(nn+1:end)), xedges);
    M = M + (null(ii)>=obs);
    % every 50 permutations, check significance
    % check requires #{null>obs} > ~10
    if ~mod(ii,50) && M >= 10
        stop = check_significance( obs, null, 0.05 );
    end
    ii = ii + 1;
end
% # perms executed
N = ii-1;

% compute p-value
if M>0
    pseudocount = 0;
else
    pseudocount = 1;
end
pval = (pseudocount + sum(null>=obs)) / N;
mdiff = median(y) - median(x);
score = obs .* sign(mdiff) .* (1 - pval );

function emd = emd_calc( x, y, xedges )
% Earth mover distance
% bin the data
x = histc( x, xedges );
y = histc( y, xedges );
% CDF
x = cumsum( x ./ sum(x) );
y = cumsum( y ./ sum(y) );
% emd
emd = sum( abs( x - y ) );

function stop = check_significance( obs, null, alpha )
% Estimate confidence bounds on p-value using normal approximation to
% binomial distribution
% Only valid if #{null>=obs} > ~10
% stop = true if 95% CI includes alpha
% For details, see Knijnenberg et al., Bioinformatics 2009
null(isnan(null)) = [];
N = length(null);
p_ecdf =  sum( null >= obs ) / N;
sigma = p_ecdf*(1-p_ecdf)/N;
lb = p_ecdf - 1.96*sigma;
% stop if bottom of CI is greater than alpha
stop = lb > alpha;
    
