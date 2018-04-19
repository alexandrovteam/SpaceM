function A = normc(A)
    % normalize per column
    maxA = max(A,[],1)
    minA = min(A,[],1)

    bsxfun(@times, bsxfun(@minus, A, minA), 1./abs(maxA - minA))
end