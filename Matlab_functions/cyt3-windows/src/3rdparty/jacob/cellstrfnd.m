function ix = cellstrfnd(cellarray, string)
% ix = cellstrfnd(cellarray, string)
% Give index as logical vector of cells in cellarray that match string or
% contain string as a substring
%@cellarray must be a vector
%@string must have the correct CASE

m = max(size(cellarray));
ix = false(1,m);

 for i = 1:m
%     %case insensitive version
    a = lower(cellarray{i});
    b = lower(string);
    if ~isempty(strfind(a,b))
%     if ~isempty(strfind(cellarray{i},string))
        ix(i) = true;
    end
end