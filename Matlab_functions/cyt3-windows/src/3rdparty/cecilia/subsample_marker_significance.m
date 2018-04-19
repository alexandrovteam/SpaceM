function subsample_marker_significance(data, current_cluster, cluster_channel, m_subsamps, channel_names)
%subsample_marker_significance(session_data(gate_context, selected_channels), current_cluster, session_data(gate_context, cluster_channel), 1000, channel_names(selected_channels))

%function for subsampling n points from gate m times and looking at L2 distances.
%n = number of cells in current cluster
%find likelihood of L2 distance for current cluster for marker

    n_cells = sum(cluster_channel == current_cluster);
    
    L2_dist = zeros(m_subsamps,size(data,2));   % m_subsamps by number of selected channels
    
    %sub sampling
    for i=1:m_subsamps,
        inds = randsample(size(data,1),n_cells);
        
        
        for j=1:size(data,2),   %looping through channels
            [P_total,IX] = ksdensity(data(:,j)); %distribution for whole population for marker
            [P_subsamp, ~] = ksdensity(data(inds,j), IX);   %distribution for subsample for marker
            
            L2_dist(i,j) = norm(P_total - P_subsamp);
        end
    end
    
    %plotting distributions
    for j=1:size(data,2),   %looping through channels
        [f,ix] = ksdensity(L2_dist(:,j));   %estimating density of L2 distributions of subsamples
        
        %figure('visible','off');
        figure();
        p = plot(ix,f);
        file_name = strcat('/Users/centh/Desktop/Columbia/Figures/from_cyt/',char(channel_names(j)),'_cluster',mat2str(current_cluster),'_subsample_L2_distribution');
        saveas(p, file_name,'tif')
    end
end