function sortedclusterheatmap( data, labels, names , rescale)
    % sortedclusterheatmap( data, labels, names )
    % data: N x M
    % labels: N x 1 or 1 x N vector of cluster assignments
    % 1 x M or M x 1 cell array of column names

data( labels == 0, :) = [];
labels( labels == 0 ) = [];
clusterSize = arrayfun(@(x) sum(labels==x), 1:max(labels) );

% sort each dimension in each cluster
mtx = [];
for c = 1:max(labels)
    
    cluster = data( labels==c, : );
    submtx = nan(size(cluster));
    for j = 1:size(cluster,2);
        submtx(:,j) = sort( cluster(:,j), 'descend' );
    end
    mtx = [mtx; submtx];
    
end

% % reorder columns for easier viewing
% mnmtx = zeros( max(labels), size(data,2));
% for c = 1:max(labels)
%     mnmtx(c,:) = mean(data(labels==c,:));
% end
% [~,co] = sortmtx( mnmtx );
% mtx = mtx(:,co);
% names = names(co);

if rescale
    mtx=mynormalize(mtx, 100);
%     mtx = bsxfun(@rdivide, mtx, prctile(mtx,99.9));
% 
%     for i=1:length(mtx)
%         ind=find(mtx(i,:)>1);
%         for j=ind(1,:)
%             mtx(i,j)=1;
%         end
%     end
end



imagesc( mtx' )
lines = cumsum( clusterSize ) + .5;
for i = 1:length(lines)-1 % don't plot last line
    hold on;
    plot( [lines(i) lines(i)], ylim(), '--', 'LineWidth', 1 ,'Color', 1*ones(1,3));
    hold off;
end
% remove empty clusters
remove = clusterSize==0;
lines(remove) = [];
k = length(unique(labels));
% set xticks in the middle of each cluster
midpoints = [];
lines = [0 lines];
for i = 1:k
    midpoints(i) = (lines(i) + lines(i+1))/2;
end
set(gca,'xtick',midpoints, 'xticklabel', unique(labels))
xlabel('cluster ID', 'fontsize', 14 )
% set(gca,'xtick',[])
if exist( 'names', 'var' )
    set(gca,'ytick',1:length(names))
    set(gca,'yticklabel', names, 'fontsize', 14)
end
colorbar;
