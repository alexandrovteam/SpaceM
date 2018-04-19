function spdists = kNNgraph( data, k, distance )

[n m] = size(data);

if ~exist('distance', 'var')
    distance = 'euclidean';
end

NS = ExhaustiveSearcher(data, 'distance', distance );

spdists = sparse( n, m );

[IDX D] = knnsearch( NS, data, 'k', k+1 );
IDX(:,1) = []; D(:,1) = []; %remove self-distances

for i = 1:n
    
    spdists( i, IDX(i,:) ) = D(i,:);
    
end

v = spdists(:);
v( v == 0 ) = [];
cutoff = full( mean(v) + 1.5*std(v) );
spdists( spdists > cutoff ) = 0;
spdists( spdists > 0 ) = 1;

