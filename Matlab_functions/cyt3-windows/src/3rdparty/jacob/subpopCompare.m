function subpopCompare( referenceData, subpopData, channelnames )
% visualize density estimates of each dimension of data, comparing
% subpopData(mxd) to referenceData(nxd), m < n.
% place subplots in order of decreasing difference as measured by histDist
% function (euclidean difference between CDFs)

[~, d] = size( referenceData );

if nargin < 3
    channelnames = cell(1,d);
end

for j = 1:d
    dists(j) = histDist( referenceData(:,j) , subpopData(:,j) );
end
[~,ix] = sort( abs(dists), 'descend' );

[nrows ncols] = numplots( d );
figure
c = 1;
for j = ix;
    subplot( nrows, ncols, c )
    dplot( referenceData(:,j) )
    dplot( subpopData(:,j) )
    title( channelnames{j} )
    box on
    set(gca,'ytick',[])
    set(gca,'xtick',0:4:8)
    xlim([-2 8])
    c = c+1;
end