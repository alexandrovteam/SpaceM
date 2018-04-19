%==============================================================================
%    mygscatter.m
%    Copyright (C) 2007  Amoolya H. Singh, singh@embl.de
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%===============================================================================
function mygscatter(X,Y,Grp,colorm, marker, leglabels, outline)
%------------------------------------------------------------------------------
%
% mygscatter(X,Y,Grp,axislabels,leglabels,outline) is an
% improvement on the Matlab program gscatter to add labels and
% outline the data points with a chosen color.  Please also see
% Matlab documentation for gscatter.
%
% Input: 
%   X, x-coordinates of data points to be plotted
%   Y, y-coordinates of data points to be plotted
%   Grp, grouping vector on data
%   axislabels, labels for X and Y axes
%   leglabels, labels for legend (text explanation of Grp)
%   outline, 1 to get light gray outline on data points, 0 o/w
%
%------------------------------------------------------------------------------

grpindices = unique(Grp)
ngroups = length(grpindices);

symbols = 'osdv^<>osdv^<>osdv^<>';
n = length(X);
h = [];
hold on
for i=1:n
  if(outline)
    c = [0.7,0.7,0.7];
    m = symbols(Grp(i));
    s = 7;
  else
    c = colorm(Grp(i),:);
    m = marker;
  end
  
  h(i) = plot(X(i),Y(i),marker,...
	      'MarkerFaceColor',colorm(Grp(i),:),...
	      'MarkerEdgeColor',c,...
	      'Marker',m);
end

ii = []; str = [];
for i=1:ngroups
  j=grpindices(i);
  ind = find(Grp==j);
  ii = [ii;ind(1)];
end

legend(h(ii), leglabels, 'Location','Best');
legend boxoff
xl = char(axislabels(1)); yl = char(axislabels(2));

hold off

end