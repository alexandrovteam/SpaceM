function [L2_density_acceptance] = L2_marker_significance(L2_dist)

%L2 dist is a matrix of L2 distances for each marker for each cluster
    
    L2_density_acceptance = zeros(size(L2_dist));
    
    for i=1:size(L2_dist,2),  %looping though channels
        
        %finding cutoff
        [f,ix] = ksdensity(L2_dist(:,i));
        
        %plotting distributions
%         figure('visible','off');
%         plot(ix,f);
%         saveas(gca, strcat(out_dir,mat2str(i),'_L2_distribution'),'tif')
        
        [~,IMAX,~,IMIN] = extrema(ksdensity(L2_dist(:,i)));

        imax = ix(IMAX(1)); %index of global maximum
        imins = ix(IMIN);   %indexes of minima
        cutoff = min(imins(imins>imax));    %cutoff is minimum immediatly to the right of global maximum
        
        L2_density_acceptance(:,i) = L2_dist(:,i) > cutoff;
      
    end
    
end