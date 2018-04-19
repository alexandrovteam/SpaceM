%% Using PhenoGraph to get meta clusters

function [cluster_mapping, centroids, meta_cluster_channel] = PhenoGraph_meta_clusters(data, cluster_channel, inds, nKNeighbors, metric)

%returns matrix where each row is a cluster and 
%col1 = meta cluster number
%col2 = cluster number
%col3 = sample/gate number
%col4 = within sample cluster number
%col5 = percentage of cells in sample belonging to cluster

    cluster_mapping = [];
    centroids = [];
    k = 1;
    cluster_sizes = [];
    
    for i=1:size(inds,2), %looping through selected gates (samples)
        cluster_mapping = vertcat(cluster_mapping, horzcat(repmat(i,[size(unique(cluster_channel(inds{i}))),1]),unique(cluster_channel(inds{i}))));
        unique_clusters = unique(cluster_channel(inds{i}));
        
        %for j=1:size(unique(cluster_channel(inds{i})),1), %looping through clusters in gate/sample
        for j=1:length(unique_clusters), %looping through clusters in gate/sample

            d = data(inds{i},:);
            %centroids(k,:) = mean(data(cluster_channel(inds{i}) == j,:));
            centroids(k,:) = mean(d(cluster_channel(inds{i}) == unique_clusters(j),:));

            %cluster_sizes(k) = size(data(cluster_channel(inds{i}) == j,:),1) / size(cluster_channel(inds{i}),1);
            cluster_sizes(k) = size(d(cluster_channel(inds{i}) == unique_clusters(j),:),1) / size(cluster_channel(inds{i}),1);

            k=k+1;
        end
    end
    
    cluster_mapping = horzcat(linspace(1,size(cluster_mapping,1),size(cluster_mapping,1))', cluster_mapping, cluster_sizes');   %col1 = cluster number, col2=sample number, col3=withing sample cluster number
  
    [meta_clusters, ~] = phenograph(centroids, nKNeighbors,'distance',lower(metric));

%     % create a sparse KNN matrix
%     sparse_adjacency_matrix = spdists_knngraph(centroids, nKNeighbors, lower(metric), 5000);  %why 5000?
% 
    
%     %finding level with highest modularity
%     mod_high = mod(end);
%     meta_clusters = cmty{find(mod == mod_high, 1,'first')};
     cluster_mapping = horzcat(meta_clusters, cluster_mapping);
    
    %creating meta_cluster_channel
    meta_cluster_channel = [];

    for i=1:size(inds,2),   %looping through gates
        met = zeros(length(inds{i}),1);
        minor_clusters = cluster_mapping(cluster_mapping(:,3) == i, 4);
        for j=1:sum(cluster_mapping(:,3) == i),
            %cluster_mapping(cluster_mapping(:,3) == i, 4),    %looping through minor clusters (can not assume that these are continous)
            met(cluster_channel(inds{i}) == minor_clusters(j)) = cluster_mapping(cluster_mapping(:,3)==i & cluster_mapping(:,4) == minor_clusters(j),1);
        end
        meta_cluster_channel = vertcat(meta_cluster_channel, met);
    end
    
end