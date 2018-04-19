function plotyy_fixZoom(AX)
% this function fixes problems with plotyy tickmarks
set(AX(2), 'XTickLabel','','XAxisLocation','Top')  %delete the xTicklabels of second x-axes and put it to top to make closed box-graph (we have to set 'Box' to 'off' in the following to remove TickMark-crosstalk between the two y-axes 
set(AX(1),'Box','off'); 
set(AX(2),'Box','off');
set(AX(1),'YTickMode','auto'); 
set(AX(2),'YTickMode','auto');
end