function [COMTY, edge_distances] = graph_louvain_cluster(adjacency_matrix, k_neighbors, varargin)
%make an adjacency matrix from this data 

    adjacency_matrix = abs(adjacency_matrix);

     [r,c] = size(adjacency_matrix);
           
            
     [Distances,I] = sort(adjacency_matrix,2,'ascend');

     %adjacency_matrix = zeros(num_events,num_events);        
     %i_sparse = zeros(1,r*k_neighbors);
     %j_sparse = zeros(1,r*k_neighbors);
     %s can just be a scalar
     i_sparse = [];
     j_sparse = [];
     edge_distances = [];       
     optargin = size(varargin,2);
     
     
     current_index=1;       
     for i=1:r
                
        for j=1:k_neighbors
            %because 1 will be itself with distance 0   
            edge_sink = I(i,j+1);
            edge_distance= adjacency_matrix(i,edge_sink);
            %adjacency_matrix(i,edge_sink) = 1; 
            if(optargin>0 && varargin{1}>0)
                if(varargin{1}<edge_distance)
                    continue;
                end
            end
            
            i_sparse(current_index) = i;
            j_sparse(current_index) = edge_sink;
            edge_distances(current_index)=edge_distance;
            current_index = current_index+1;
            
        end
     end

     sparse_adjacency_matrix = sparse(i_sparse,j_sparse,1, r, r);
     
    [COMTY ending] = cluster_jl(sparse_adjacency_matrix,1,1,1,1);
    
end

