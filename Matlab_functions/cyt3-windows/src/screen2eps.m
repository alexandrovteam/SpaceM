% SCREEN2EPS Generate an EPS file of the current figure with
%   dimensions consistent with the figure's screen dimensions.
function screen2eps(h,filename)

% create directory if needed
[pathstr, ~, ~] = fileparts(filename);
if numel(pathstr) >0 && ~exist(pathstr, 'dir')
   mkdir(pathstr); 
end

% print
print(h,'-dpdf', filename, '-r0');

end
