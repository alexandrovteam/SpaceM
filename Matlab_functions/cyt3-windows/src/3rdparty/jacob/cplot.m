function [rho percentnz] = cplot(x1, x2, pos)
%[rho percentnz] = cplot(x1,x2,pos)
%generate contour plot for two-dimensional data
%return spearman's RHO for data and percent of data for which all(data) > 0
%option to exclude negative pairs of values
%@PARAMS
% x1: mx1 or mx2 matrix
% x2: mx1 vector or empty if x1 is mx2
% pos: if pos = 'pos', negative pairs of values are excluded

%nnv results in plotting only the pairs of valus without non-negative
%values. enter 1 for this option

inputsize = size(x1,1);

if nargin > 2
    pos = true;
else
    pos = false;
end

if size(x1,2) == 2
    if ~exist('x2','var')
        x2 = 0;
    end
    pos = x2;
    x2 = x1(:,2);
    x1 = x1(:,1);
end

% if size(x1,2) == 2
%     if ~exist('x2','var')
%         x2 = 0;
%     end
%     nnv = x2;
%     x2 = x1(:,2);
%     x1 = x1(:,1);
% end
%
% if strmatch(nnv, 'nnv')
%     idx1 = x1 > 0;
%     idx2 = x2 > 0;
%     x1 = x1(idx1 & idx2);
%     x2 = x2(idx1 & idx2);
% end

%if ischar(pos) && strmatch(pos,'pos')
if pos
    [~,ix] = pospairs([x1 x2]);
    x1 = x1(ix);
    x2 = x2(ix);
else
    percentnz = NaN;
end

[~, density, x, y] = kde2d([x1 x2], 256);
%figure
contour(x, y, density, 12);

%rho = corr(x1, x2, 'type', 'spearman');

if ~exist('percentnz','var')
    percentnz = numel(x1)/inputsize;
end
end