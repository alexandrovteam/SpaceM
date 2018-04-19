% returns true if an input string ends with the given endstr string; false
% otherwise.
function ew=endswith(str, endstr)

inds = strfind(str, endstr);

% index found and it corresponds with the end of the string
ew= ~isempty(inds) && inds(end) == numel(str) - numel(endstr) + 1;

end