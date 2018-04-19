function h = scatterlabels( map, labels, markersize )
% scatterlabels( map, labels, markersize )
n = length(labels);
if n ~= length(map)
    error('Map and labels must have same number of rows')
end
% variable marker size
if nargin < 3
    markersize = 10*ones(n,1);
elseif length(markersize) == 1
    markersize = markersize .* ones(n,1);
end

% figure
hold all

% number of groups
ul = unique(labels);
numgroups = length(ul);
ul = reshape(ul,1,numgroups);
cmap = distinguishable_colors(numgroups);

% are labels numeric categories?
numericlabels = isnumeric(labels);
h = [];
if numericlabels
    
    % if one group is 0, change colors
    if any(ul==0)
        cmap = distinguishable_colors(numgroups-1);
        cmap = [.5 .5 .5; cmap];
    end
    
    L = cell(1,numgroups);
    tick = 1;
    for i = ul
        ix = labels == i;
        h(end+1) = scatter(map(ix,1),map(ix,2),markersize(ix),cmap(tick,:),'fill');
        L{tick} = num2str(i);
        tick=tick+1;
    end
    legend( L, 'location', 'northeastoutside' )
    
else
    
    for i = 1:numgroups
        ix = strcmp( labels, ul{i} );
        h(end+1) = scatter(map(ix,1),map(ix,2),markersize(ix),cmap(i,:),'fill');
    end
    legend( ul, 'location', 'northeastoutside' )
    
end
        