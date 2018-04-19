% Function for finding L2 distances between marker distribution for cluster
% and marker distribution for all clusters for a single gate. Sign (-/+)
% then determined based on whether the median for the cluster for the
% marker is higher or lower than the median for that marker across all
% clusters

function [hl_L2_dist] = L2_dist_heatmap_data(gate_data, cluster_channel)

    %testing if the number of clusters is reasonable
    if length(unique(cluster_channel)) > 500,
        error('Number of clusters is too high to run L2_dist_heatmap_data\n');
        
    end
        
    gate_medians = median(gate_data);   %finding medians for each marker for cells belonging to gate

    n_clusters = length(unique(cluster_channel));  %finding number of clusters in gate
    hl_median = zeros(n_clusters, length(gate_medians)); %defining matrix for storing above/below median information (NxM, N=clusters, M=markers)

    clusters_in_sample = unique(cluster_channel);

    %determening if median for cluster for for each channel is higher/lower than median across clusters
    for i=1:length(unique(clusters_in_sample)),  %lopping through clusters
        for j=1:size(gate_data,2),     %looping through markers
            m = median(gate_data(cluster_channel == clusters_in_sample(i), j)); %median for marker in cluster
            if m < gate_medians(j),
                t = -1;
            elseif m == gate_medians(j),
                t = 0;
            elseif m > gate_medians(j),
                t = 1;
            end
            hl_median(i,j) = t;
        end
    end


    % L2 as measure of distance between cluster and whole population
    L2_dist = zeros(n_clusters, size(gate_data,2));
    
    for j=1:size(gate_data,2),  %looping through each marker

        %distribution for whole population for marker
        [P_total,IX]=ksdensity(gate_data(:,j));

        %loop through each cluster and find distance to whole population
        for i=1:length(clusters_in_sample),
            [P_cluster]=ksdensity(gate_data(cluster_channel==clusters_in_sample(i),j),IX);
            L2_dist(i,j)=norm(P_total-P_cluster, 2);
        end
    end


    %L2 combined with higher/lower than median
    hl_L2_dist = L2_dist .* hl_median;
      
end
