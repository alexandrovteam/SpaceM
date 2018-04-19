function density = kde2d_color_hist(data, varargin)
%KDE2D_COLOR_HIST Draws a kerenel density plot along y axes.
%   input args:
%   data            - NX2 data matrix
%   'thresh'        - a minimum under which not to normalize (default=0.0001)
%   'axis'          - normalize along x (2 - along 2nd dim) or y (1 - along 1st dim) axis. specify the dim.
% 
%   output
%   density         - the density matrix normalized along y axes

if (size(data,2) ~= 2)
    error('data has to be an N by 2 array where each row represents a two dimensional observation')
end

norm_dim_thresh = 0.0001;
axis_norm = 0;
for i=1:length(varargin)
    switch lower(varargin{i})
        case 'axis'
            axis_norm = varargin{i+1};
            if ~(axis_norm==1 || axis_norm == 2)
                error('normalize along dim 1 or 2');
            end
        case 'thresh'
            norm_dim_thresh = varargin{i+1};
    end
end

hist_l = 10;
dev_l = 1;

% get density matrix
[~, density, x, y] = kde2d(data, 128);

% normalize density matrix (?? why isn't it returned normalized ??)
density = mat2gray(density, [0 max(density(:))]);

max_dens_y = max(density, [], 1);
max_dens_x = max(density, [], 2);
max_dens = [max_dens_y(:) max_dens_x(:)];

if axis_norm
    % ignore maximum values that are to low for numerical reasons
    max_dens_for_norm = max_dens(:,axis_norm);
    max_dens_for_norm( max_dens(:,axis_norm) < norm_dim_thresh) = inf;

    % devide each column by the maximum of that column (normalize col-wise)
    %%%% TODO !! density_bsx = bsxfun(@rdivide,density,max_dens_for_norm);
    
    %%% temp, ineffiect way but works.
    if axis_norm == 2
        density = density';
    end
    for i=1:size(density,1) % density is an nXn matrix
        
        density(:, i) = density(:, i)/max_dens_for_norm(i); %cols (y)
    end
    
    if axis_norm == 2
        density = density';
    end
end
    
% add 1d histogram\intesity levels along x and y axis
density = [repmat(max_dens_x(:), 1, hist_l) density];
density = [repmat([nan(1, hist_l) max_dens_y(:)'], hist_l, 1); density];

% add some devider between the 1d histograms and the image
density(hist_l:(hist_l-dev_l+1), :) = NaN;
density(:, hist_l-dev_l+1:hist_l)   = NaN;

% calculate the axis ticks
xmin = min(data(:, 1));
ymin = min(data(:, 2));
xmax = max(data(:, 1));
ymax = max(data(:, 2));

added_hists = (xmax-xmin)*(hist_l/(size(density, 1)-hist_l));
xbins = xmin-added_hists:0.5:xmax;
ybins = ymin-added_hists:0.5:ymax;

% plot the hist
cla;
cmap = gray(20);
cmap = cmap(end:-1:1, :);
cmap(1:1, :) = 1;
colormap(cmap);
legend('off');
colorbar('delete');

nanimagesc(density, cmap, 'nancolor', [0.15 0.15 0.15], 'xbins', xbins, 'ybins', ybins);

colorbar; 
set(gca,'YDir','normal');

end

