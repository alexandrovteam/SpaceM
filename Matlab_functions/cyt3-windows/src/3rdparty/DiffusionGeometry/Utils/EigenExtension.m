function [EigenVecs_ext,EigenVals_ext] = EigenExtension( X, Y, G, opts )

% function EigenExtension( G, opts )
%
% IN:
%   X       : D by N set of points on which G was constructed
%   Y       : D by M set of points at which to extend the eigenvectors
%   G       : graph and diffusion information as constructed by GraphDiffusion
%   opts    : structure of options, including the following fields:
%               [kEigenVecs] : how many eigenvectors to extend, up to the number that was already computed. 
%                               Default: all the ones in G.EigenVecs
%
% OUT:
%   EigenVecs_ext   : (N+M) by opts.kEigenVecs extended eigenvectors
%   EigenVals_ext   : opts.kEigenVecs 
%

if nargin<4, opts = []; end;
if ~isfield(opts,'kEigenVecs'), opts.kEigenVecs = size(G.EigenVecs,2); else
                                opts.kEigenVecs = max([min([opts.kEigenVecs,size(G.EigenVecs,2)]),1]);end;

                            
%% Add points to graph
% Find nearest neighbors of new points among the old points
[count,idxs,dists,NNInfo] = nrsearch( X, Y, G.Opts.kNN, 0, struct('ReturnAsArrays',true,'NNInfo',G.DistInfo.NNInfo) );
i_new = (1:length(Y))';i_new=i_new*ones(1,size(idxs,2));i_new=i_new';i_new=i_new(:);i_new=i_new;
j_new = idxs';j_new=j_new(:);
for i = 1:size(dists,1),
    dists(i,:) = dists(i,:)/dists(min([size(dists,2),G.Opts.kNNAutotune]));
end;
dists = dists'; dists=dists(:);
T_new = sparse(i_new,j_new,dists,size(Y,2),size(Y,2));

% Find nearest neighbors of new points among the new points
[count2,idxs2,dists2,NNInfo2] = nrsearch( Y, Y, G.Opts.kNN, 0, struct('ReturnAsArrays',true) );
i_new2 = (1:length(Y))';i_new2=i_new2*ones(1,size(idxs2,2));i_new2=i_new2';i_new2=i_new2(:);i_new2=i_new2;
j_new2 = idxs2';j_new2=j_new2(:);
for i = 1:size(dists2,1),
    dists2(i,:) = dists2(i,:)/dists2(min([size(dists2,2),G.Opts.kNNAutotune]));
end;
for i = 1:length(i_new),
    s_new(2i) = exp(-dists2(i_new2(i),j_new2(i))*dists(i_new2(i),j_new2(i)));
end;


T_new2 = sparse(double(i_new2),double(j_new2),s_new2,size(Y,2),size(Y,2));
T_new2 = 0.5*(T_new2+T_new2');

npm = size(G.T,1)+size(Y,2);

[i,j,s] = find(G.T);

T = sparse( double([i;i_new+size(G.T,1);j_new]),double([j;j_new;i_new+size(G.T,1)]),double([s;s_new;s_new]),npm,npm );
                            
% Extend the eigenvector to enlarged graphs                            
[EigenVecs_ext,EigenVals_ext] = eigs( T, opts.kEigenVecs, 'LM', struct('disp', 0,'isreal',true,    ...
                                        'p',2*opts.kEigenVecs, ...
                                        'v0',[G.EigenVecs,randn(size(G.EigenVecs,1),opts.kEigenVecs-size(G.EigenVecs,2))]));

EigenVals_ext = diag( EigenVals_ext );

return;