function scatter_discrete(x, y, gcolors, dot_size)
        
    groups = unique(gcolors);
    nGroups = numel(groups);
    
    %color_per_group = distinguishable_colors(nGroups);
    
    %only 32 colors
    map=[[0,0,1];            [0.929,0.694,0.125];[0.166,0.374,0.088];
         [0.482,0.408,0.933];[0.545,0.270,0.074];[0.498,1.000,0.831];
         [0,0.545,0.545];    [0,0.247,0.541];    [0.604,0.804,0.196];
         [0.867,0.627,0.800];[1,0,1];            [0.494,0.184,0.556];
         [1.000,0.078,0.576];[1.000,0.855,0.725];[0.561,0.737,0.561];
         [1.000,0.549,0];    [0.502,0.502,0];    [0,0.447,0.741];
         [0.85,0.325,0.098]; [0.466,0.674,0.188];[0.201,0.845,0.933];
         [1,0,0];            [0,1,0];            [0.194,0.684,0.256];
         [0.65,0.125,0.098]; [0.729,0.494,0.025];[1,1,0];            
         [0.2,0.647,0.941];  [0.929,0.894,0.525];[0.282,0.239,0.545];
         [0.766,0.974,0.388];[0.201,0.945,0.533]];

    color_per_group = colormap(map);
    
    % plot each color in a seperate group (otherwise the legend won't work)
    for i=1:nGroups
        curr_group = groups(i);
        curr_color_inds = (curr_group==gcolors);

        scatter(x(curr_color_inds),y(curr_color_inds),...
                dot_size(curr_color_inds),...
                color_per_group(i, :), 'fill');

        hold on;      
    end
    
	% find axis limits
    x_range  = max(x)-min(x);
    x_buffer = x_range*.03;
    x_lim    = [min(x)-x_buffer, max(x)+x_buffer];
    
    y_range  = max(y)-min(y);
    y_buffer = y_range*.03;
    y_lim    = [min(y)-y_buffer, max(y)+y_buffer];

    xlim(x_lim);
    ylim(y_lim);

end