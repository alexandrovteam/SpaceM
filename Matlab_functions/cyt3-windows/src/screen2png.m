% SCREEN2JPEG Generate a JPEG file of the current figure with
%   dimensions consistent with the figure's screen dimensions.
function screen2png(h,filename)

% create directory if needed
[pathstr, ~, ~] = fileparts(filename);
if numel(pathstr) >0 && ~exist(pathstr, 'dir')
   mkdir(pathstr); 
end

set(h,'Units','pixels');

pos = get(h,'Position');
newpos = pos/100;

% sets the position\size of the current graphical object before printing
set(h,'PaperUnits','inches',...
     'PaperPosition',newpos)

% print
print(h,'-dpng', filename, '-r200');
drawnow;
end
