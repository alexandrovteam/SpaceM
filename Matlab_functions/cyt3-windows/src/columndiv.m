% divie each column by the value in vec in the same column position
function mat = columndiv(mat, vec)
    mat = bsxfun(@rdivide,mat, vec);
end