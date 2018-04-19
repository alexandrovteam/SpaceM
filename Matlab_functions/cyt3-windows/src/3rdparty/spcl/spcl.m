function [clusters, evalues, evectors] = spcl(data, nbclusters, varargin)
%
%spcl is a spectral clustering function to assemble random unknown data
%into clusters. after specifying the data and the number of clusters at
%first, next parameters are how to construct the similarity graph (the
%default is the fully connected 'fully'), what laplacian matrix to construct (default
%is the 'unormalized', what eigen vectors to use to sort data and what
%algorithm to sort with (default is eigen vector '2' and 'pn' sorting).
%
%data is the unknown data, its a matrix where each row represent a data
%vector and its length is the number of samples
%
%nbclusters is the number of clusters used to assemble data
%
%similarity graph will be constructed as fully connected graph
%first parameter is the name of the function to use, the second is
%the parameter to pass to the function (if the parameter function name
%is not found the default will be gaussian function and the parameter found
%will be its sigma)
%
%laplacian matrix parameters:
%'unormalized' - unormalized laplacian matrix
%'sym' - normalized symetric laplacian matrix
%'rw' - normalized assymetric laplacian matrix
%
%algorithm for organizing eigen vectors:
%'np' - one eigen vector must be used, if will put positive values in class
%1 and negative values in class 2
%'kmean' - a kmean algorithm will be used to cluster the given eigen
%vectors
%
%finally an eigen vector choice can be added, it can be one value or a
%vector [vmin vmax] or a matrix defining several intervals. if not found
%the default will be 2

plotchoices = {'bo','r+','md','k*','wv'};
lapmatrixchoices = {'unormalized', 'sym', 'rw'};
algochoices = {'np', 'kmean'};
func = 'gaussfunc';
count = 1;

%%get all the parameters%%%
if(ischar(varargin{count}))
        
    func = varargin{count};
    count = count + 1;
end

params = varargin{count};
count = count + 1;

if(length(varargin) >= count)
    
    if(sum(strcmp(varargin{count}, lapmatrixchoices)) == 0)

        lapmatrixchoice = 'unormalized';
    else

        lapmatrixchoice = varargin{count};
        count = count + 1;
    end

    if(length(varargin) >= count)
        
        if(sum(strcmp(varargin{count}, algochoices)) == 0)

            clusteralgo = 'np';
        else
            clusteralgo = varargin{count};
            count = count + 1;
        end

        if(length(varargin) >= count)

            eigv = varargin{count};
        else
            
            eigv = [2 2];
        end
    else
        clusteralgo = 'np';
        eigv = [2 2];
    end
else
    
    lapmatrixchoice = 'unormalized';
    clusteralgo = 'np';
    eigv = [2 2];
end
%%all parameters are got%%%
sprintf('graph choice is fully connected\nLaplacian choice is %s\nCluster algorithm is %s', lapmatrixchoice, clusteralgo)
[nbsamples, dim] = size(data);
wmat = zeros(nbsamples);

for i = 1: nbsamples - 1
    
    wmat(i, i + 1: end) = feval(func, repmat(data(i, :), nbsamples - i, 1), data(i + 1: end,:), params);
end

wmat = wmat + wmat';
dmat = diag(sum(wmat, 2));

if(strcmp(lapmatrixchoice, 'unormalized'))
    
    laplacian = dmat - wmat;
else
    if(strcmp(lapmatrixchoice, 'sym'))
        
        laplacian = eye(nbsamples) - (dmat^-0.5) * wmat * (dmat^-0.5);
    else
        if(strcmp(lapmatrixchoice, 'rw'))
            
            laplacian = eye(nbsamples) - (dmat^-1) * wmat;
        end
    end
end

[evectors, evalues] = eig(laplacian);

newspace = evectors(:, eigv(1,1): eigv(1,2));
n = size(eigv);
for i = 2: n(1)
    
    newspace = [newspace evectors(:, eigv(i,1): eigv(i,2))];
end

if(strcmp(clusteralgo, 'kmean'))
    
    clusters = kmeans(newspace, nbclusters);
else
    clusters = 1 + (newspace > 0);
end

if(dim == 2)
    figure;
    for i = 1: nbclusters

        points = data(clusters == i, :);
        plot(points(:,1), points(:,2), plotchoices{i});
        hold on;
    end
    title('clustered data using spectral clustering');
    grid on;
end