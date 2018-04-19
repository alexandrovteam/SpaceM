function myscatter_by_point(x, y, color_var, dot_size, discrete)

hold off;
    
    x_range = max(x)-min(x);
    y_range = max(y)-min(y);
    
    x_lim = [min(x)-x_range*.02,max(x)+x_range*.02];
    y_lim = [min(y)-x_range*.02, max(y)+y_range*.02];
    
    if discrete
        colors = distinguishable_colors(max(color_var));
    else
        cmapsize = 20;
        colors = jet(cmapsize);
        color_var = ceil(mynormalize(color_var, 99.5)*cmapsize);
        color_var(color_var>20) = cmapsize;
        color_var(color_var<1) = 1;
    end
    
    for i=1:size(x),

        scatter(x(i),y(i), dot_size(i), colors(color_var(i),:), 'fill');
        hold on;
        
    end
 
	xlim(x_lim);
    ylim(y_lim);

    %xlabel('tSNE1');
    %ylabel('tSNE2');
    
    if discrete
        %adding legend
        legend_handle = legend(gca,num2str(unique(color_var)), 'location', 'eastoutside');
        legend_markers = findobj(get(legend_handle, 'Children'), 'marker', 'o');    %find handles of markers in legend

        p = length(legend_markers);
        for i=1:length(legend_markers), %set color of markers in legend
            set(legend_markers(p), 'markerfacecolor', colors(i,:));
            p=p-1;
        end
    else
        colorbar
    end
    hold off;

end