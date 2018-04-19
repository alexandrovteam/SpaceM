function [L2_density_acceptance] = percentile_marker_significance(L2_dist, percentile)

%L2_dist is a matrix of L2 distances for each marker for each cluster. Percentile is the percentile cutoff wanted
    
    L2_density_acceptance = zeros(size(L2_dist));
    
    for i=1:size(L2_dist,2),  %looping though channels
        
        cutoff = quantile(L2_dist(:,i), percentile);
        L2_density_acceptance(:,i) = L2_dist(:,i) > cutoff;
      
    end
    
end