function scatter_by_point(x, y, gcolors, dot_size)
    
    % plot each color in a seperate group (otherwise the legend won't work)
    
    % get the different colors (groups)
    groups = unique(gcolors);
    nGroups = numel(groups);
    
    % create a color map
    map = distinguishable_colors(nGroups);     
%     colormap(map); % removing bc the colors are individually spec
    
    % For each group
    for i=1:nGroups
        curr_group = groups(i);
        
        % get indices of current color group
        curr_color_inds = (curr_group==gcolors);

        % scatter the current group
        scatter(x(curr_color_inds),y(curr_color_inds),...
                dot_size(curr_color_inds),...
                map(i, :), 'fill');

         hold on; 
    end
end