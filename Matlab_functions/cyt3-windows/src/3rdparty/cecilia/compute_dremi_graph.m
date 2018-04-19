
%function [m_cdata] = pairwise_marker_distance(all_cdata, n_ct_clusters, n_samples, all_clusters, time_ct, marker_channel_names, out_dir, marker_channels)
function [g, ret] = compute_dremi_graph(data, current_cluster, channel_names)
    
    channel_names = lower(channel_names);

    %min and max values for markers at time point
    miny_val = min(min(data));
    maxy_val = max(max(data));

    if size(data,1) < 300,  %could be a more intelligent cut off
        fprintf (strcat('Not enough data points to compute DREMI for cluster\t',mat2str(current_cluster),'\n'))
    else
        
        %making cdata object
        cdata = cytof_data('empty');
        cdata = cdata.add_data_matrix_header(data, channel_names);
        
        %computing DREMI
        [cdata, mi_matrix ] = cdata.pairwise_mi_compute(channel_names, .75,'MinMaxY', miny_val, maxy_val); %computing pairwise DREMI
        [cdata, sig_edges, sig_mis, sig_edge_matrix] = cdata.write_mi_graph(channel_names, mi_matrix, .25, strcat('tmp.txt'), 'choose');    %pruning
        cdata = cdata.prune_transitive_edges();
        [cdata, ret] = cdata.draw_denovo_DREMI_graph();   %drawing graph

        if ret == 1,    %if draw_denovo_DRMI_graph returns 1 (meaning that edges where found)
            dolayout(cdata.gObj);
            set(cdata.gObj,'ShowArrows', 'off'); 
            view(cdata.gObj);
            %g = biograph.bggui(cdata.gObj);
        
        else
            fprintf(strcat('No strong edges in graph for cluster\t',mat2str(current_cluster),'\n'));
        end
    end

end