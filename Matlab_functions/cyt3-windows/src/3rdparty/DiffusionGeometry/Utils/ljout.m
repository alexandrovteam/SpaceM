function ljout(string, len)
% function ljout(string, len)
%
% Print a string with left justified formatting and a given length.

if length(string) > len
   string = string(1:len);
else
   string = [string ' ' * ones(1, len-length(string))];
end

fprintf('%s', string);

return;