% PhenoGraph: Robust graph-based phenotype discovery
% -----------------------------------------------------------------------
%  [labels,communities,G] = phenograph( data, k, varargin )
%
% INPUT:
%     data - n x d matrix of n single cells
%     k    - number of neighbors for graph construction
% Options:
%     'distance', distance - string input for knnsearch
%     'G', graph
%     'IDX', IDX
%     'graphtype', 'jaccard' or 'knn'
% OUTPUT:
%     labels - n x 1 vector of cluster assignments
%     G      - n x n sparse matrix representation of the (lower triangle)
%              of the Jaccard graph used for community detection
% -----------------------------------------------------------------------
function [labels,communities,G,uniqueID] = phenograph( data, k, varargin )

fprintf(1,'Starting Phenograph on %i points in %i dimensions...\n',...
        size(data,1), size(data,2));
    
distance = 'euclidean';
graphtype = 'jaccard';
if nargin >2
    for i = 1:2:length(varargin)
        switch varargin{i}
            case 'G'
                G = varargin{i+1};
            case 'IDX'
                IDX = varargin{i+1};
            case 'distance'
                distance = varargin{i+1};
            case 'graphtype'
                graphtype = varargin{i+1};
        end
    end
end
Jaccard = strcmpi( 'jaccard', graphtype );

if ~exist('IDX','var') && ~exist('G','var')
    % Nearest neighbors
    fprintf(1,'Finding %i nearest neighbors...\n',k);
    [IDX,D] = knnsearch(data,data,'k',k+1,'distance',distance);
    IDX(:,1) = []; D(:,1) = [];
end

if ~exist('G','var') && Jaccard
    % Jaccard graph
    fprintf(1,'Building Jaccard Graph...\n');
    G = knn2jaccard(IDX);
end

if ~exist('G','var') && ~Jaccard
    % knn graph
    fprintf(1,'Building knn graph...\n');
    G = idx2knn(IDX,D);
end

if any(isinf(nonzeros(G)))
    [~,j] = find(isinf(G));
    display(unique(j));
    error('Graph contains infinite weights. Check your data for identical points');
end

work_dir = pwd;

[curr_path, ~, ~] = fileparts(mfilename('fullpath'));
curr_path = [curr_path filesep];

cd(curr_path);

try 
    fprintf(1, 'Writing graph G...\n');
    
    % Write graph to file
    Graph2Binary(G,'G');

    % Run Louvain on file for multiple iterations
    niter = 20;
    [c,Q,labels,communities] = LouvainfromBin('G.bin',niter);

    llim = max([ceil(length(labels)./1e4) 1]);
    labels = sortlabels(labels,llim);

    %generate unique id 
    uniqueID = genarateID;
        
    cd(work_dir);
catch ME
    cd(work_dir);
    rethrow(ME);
end
end

function c = sortlabels( c, cutoff )
% c = sortlabels( c, cutoff )
% Rearrange classlabels in c so the first class is the largest, etc.
% if c is a cell array, perform iteratively on each cell
% if cutoff is provided, entries in c belonging to classes with size <
% cutoff are assigned to class 0.

if nargin < 2
    cutoff = 0;
end

if iscell( c )
    for idx = 1:length(c)
        c{idx} = sortlabs( c{idx}, cutoff );
    end
    
else % c is a single unsorted labels vector
    
    c = sortlabs( c, cutoff );
    
end
end
function sl = sortlabs( labs, cutoff )

i = 1;
classSize = arrayfun(@(x) sum(labs==x), 1:max(labs));
sl = zeros(size(labs));
remove = sort( classSize <= cutoff, 'ascend');
[~,ix] = sort(classSize,'descend');
ix(remove) = [];
for assignment = ix
    sl( labs == assignment ) = i;
    i = i+1;
end
end