function cost= compare_maps(map1,map2, measureType)
% This function is using t-SNE's cumputing of matrix P and matrix Q in a 
% low-dimensional maps.

    P= compute_Pvalues(map1);
    Q= compute_Pvalues(map2);
    
    if (measureType == 1)
        cost=kldivergences(P,Q);
    else
        cost = JSdivergences(P,Q);
    end
end


function value = compute_Pvalues(map) 
% Computing joint probability that point i and j are neighbors
    n = size(map, 1);                                                      % Initialize number of instances
    sum_map = sum(map .^ 2, 2);                                            % make sure the value sum to one
    num = 1 ./ (1 + bsxfun(@plus, sum_map, ...
        bsxfun(@plus, sum_map', -2 * (map * map'))));                      % Student-t distribution
    num(1:n+1:end) = 0;                                                    % set diagonal to zero
    value = max(num ./ sum(num(:)), realmin); 
end

function cost = kldivergences(P,Q)                                         %Kullback-Leibler divergence of the maps
    const = sum(P(:) .* log(P(:)));                                        % constant in KL divergence
    cost = const - sum(P(:) .* log(Q(:)));   
    
end

function cost = JSdivergences(P,Q)                                         %Jensen-Shannon divergence of the maps
    M = 0.5*(P + Q);
    cost = 0.5*(kldivergences(P,M) + kldivergences(Q,M));   
end