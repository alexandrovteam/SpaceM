function distance = histDist(referenceData, sampleData, metric)
%distance = histDist(referenceData, sampleData, metric, P)
%default if no metric specified: euclidean

if nargin < 3
    metric = 'euclidean';
end
%set grid
m = min( [min(referenceData) min(sampleData)] );
M = max( [max(referenceData) max(sampleData)] );
interval = (M-m) / 100; %value?
xi = m:interval:M;
%compute CDFs and distance
cdfRef = ksdensity( referenceData, xi, 'function', 'cdf' );
cdfSamp = ksdensity( sampleData, xi, 'function', 'cdf' );

distance = pdist2(cdfRef, cdfSamp, metric);


distance = distance * sign( median(sampleData) - median(referenceData) );

end