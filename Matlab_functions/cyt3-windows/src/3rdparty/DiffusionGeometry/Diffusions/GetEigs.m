function [V, D] = GetEigs(T, k, P, Opts)

% function [V, D] = GetEigs(T, k, P)
%
% GETEIGS fetches the k largest magnitude eigenvales for the matrix T.
% [Opts] : structure that contains the following parameters:
%           'TakeDiagEigenVals' : if 1, returns the eigenvalues as a vector rather than as a diagonal matrix. Default: 0.
%
%

if nargin<4,
    Opts = [];
end;
if ~isfield(Opts,'TakeDiagEigenVals'),
    Opts.TakeDiagEigenVals = 0;
end;
   

if ~iscell(T),
    if (~exist('P')) | (isempty(P)),
        P = speye(size(T));
    end

    [V, D] = eigs(T, max([1,min([k,size(T,1)-1])]), 'LM', struct('disp', 0,'tol',1e-4,'maxit',1000));

    % sort them
    v = diag(D);
    [v idxs] = sort(v, 'descend');

    D = D(idxs, idxs);
    V = V(:, idxs);
    V = P*V;
    for k = 1:size(V,2),
        V(:,k) = V(:,k)/norm(V(:,k));
    end;
    
    if Opts.TakeDiagEigenVals,
        D = diag(D);
    end;
    
    %norm(full(T*V - V*D))
else
    for lk = 1:length(T),
        if (~exist('P')) | (isempty(P{lk})),
            P{lk} = speye(size(T{lk}));
        end

        [lV, lD] = eigs(T{lk}, max([0,min([k,size(T,1)-1])]), 'LM', struct('disp', 0));

        % sort them
        v = diag(lD);
        [v idxs] = sort(v, 'descend');

        if Opts.TakeDiagEigenVals,
            D{lk} = diag(lD(idxs,idxs));
        else
            D{lk} = lD(idxs, idxs);
        end;
        V{lk} = lV(:, idxs);
        V{lk} = P{lk}*V{lk};
        V{lk} = V{lk}/norm(V{lk});
    end;
end;

