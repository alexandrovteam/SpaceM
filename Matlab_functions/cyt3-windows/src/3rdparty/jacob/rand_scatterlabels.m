function h = rand_scatterlabels( map, labels, markersize, cmap )
% scatter enties in Nx2 map with colored labels and markersize
% randomize plotting order to improve visualization

n = length(labels);
if n ~= size(map,1)
    error( 'map and labels must have same length' );
end

% variable marker size
if nargin < 3
    markersize = 10*ones(n,1);
elseif length(markersize) == 1
    markersize = markersize .* ones(n,1);
end

if nargin < 4
    gencmap = true;
else
    gencmap = false;
end

hold all

% number of groups
ul = unique(labels);
numgroups = length(ul);
ul = reshape(ul,1,numgroups);
if gencmap
    cmap = distinguishable_colors(numgroups);
end
    
C = zeros(n,3);
% for g = 1:numgroups
q = 1;
for g = ul
    ix = labels==g;
    C(ix,:) = repmat(cmap(q,:),sum(ix),1);
    q = q+1;
end

rix = randperm(n);
h = zeros(1,numgroups);
for i = rix
    g = labels(i);
    idx = find(ul==g);
    h(idx) = scatter(map(i,1),map(i,2),markersize(i),C(i,:),'fill');
end
