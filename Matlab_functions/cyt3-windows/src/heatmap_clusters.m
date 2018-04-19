function heatmap_clusters(data, labels, label_names)

[~,ix] = sort(labels);
data( data < 0 ) = 0;
imagesc( data(ix,:)' )
colorbar
set(gca,'ytick',1:size(data, 2))
set(gca,'yticklabel',label_names)

end