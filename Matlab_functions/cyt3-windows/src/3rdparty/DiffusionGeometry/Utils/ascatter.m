function h = ascatter(points, f)
% function ascatter(points, [f])
%
% ASCATTER is a frontend for the MATLAB scatter and scatter3 commands.
% Given a list of points n either R^2 or R^3, it will produce a scatter
% plot.  If the optional argument f is specified, the values of f at
% each point will be used as a color value, producing a representation
% of the function f.
%
% In:
%    points = an dxn matrix specifying n points in R^d
%    f      = an optional argument
%
% Out:
%    h      = handle for the scatter plot figure
%
% Dependencies:
%    none
%
% Version History:
%    jcb       5/2005         created by fixed up plotplot.m
%    jcb       8/2006         points matrix transposed

[D,N] = size(points);


if D==2
   if exist('f')
      h = scatter(points(1,:), points(2,:), repmat(200, N, 1), f, '.');
   else
      h =  scatter(points(1,:), points(2,:), repmat(200, N, 1), '.');
   end
elseif D==3
   if exist('f')
      h = scatter3(points(1,:), points(2,:), points(3,:), repmat(500, 1, N), f, '.');
   else
      h = scatter3(points(1,:), points(2,:), points(3,:),repmat(500, 1, N),'.');
   end
else
   fprintf('ascatter.m: Points must be in R^2 or R^3\n');
   h = [];
end