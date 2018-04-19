function SPR_distances = SPR_dist_heatmap(data, cluster_channel, niter)

    clusters_in_sample = unique(cluster_channel);
    SPR_distances = zeros(length(clusters_in_sample), size(data,2));

    %looping through each marker
    for i=1:size(data,2)

        pop_marker_dist = data(:,i);

        %looping through each cluster        
        for j=1:length(clusters_in_sample)

            cluster_marker_dist = data(cluster_channel == clusters_in_sample(j),i);

            [score,pval,emd,mdiff,null] = SARA(pop_marker_dist, cluster_marker_dist, niter);
            
            SPR_distances(j,i) = score;
        end
    end
end