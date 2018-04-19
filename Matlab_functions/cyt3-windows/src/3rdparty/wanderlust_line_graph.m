function [ pts, pt_indices, pt_mean_intensity, pt_std_intensity ] = wanderlust_line_graph( wanderlust, data, legend_labels, num_pts, window_width, v )
% [ pts, pt_indices, pt_mean_intensity, pt_std_intensity ] = wanderlust_line_graph( wanderlust, data, markers, num_pts, window_width )
%
% plot a line graph of marker intensities across the trajectory.
%
% (for each of num_pts equidistant points across wanderlust, pick all points in window_width
% diameter, use these to plot)
%  
% TODO: rewrite this travesty of a comment.
group_ids = unique(v);
if (length(group_ids) > 1) 
    legend_labels = remove_repeating_strings(legend_labels);
    legend_labels{end+1} = 'average';
end

% normalize each marker to 0-1 range using top and bottom one percentiles
% norm_bottom = prctile( data, 2 );
% norm_top = prctile( data, 98 );
% data = ( data - repmat( norm_bottom, length( data ), 1 ) ) ./ ( repmat( norm_top, length( data ), 1 ) - repmat( norm_bottom, length( data ), 1 ) );
% data = max(data, 0);
% data = min(data, 1);

%normalize the wanderlust channel?? but whyyyyhhhhaaahhhyyyy?
norm_bottom = prctile( wanderlust, 2 );
norm_top = prctile( wanderlust, 98 );
wanderlust(wanderlust<norm_bottom) = norm_bottom;
wanderlust(wanderlust>norm_top) = norm_top;
wanderlust = (wanderlust-min(wanderlust));
wanderlust = wanderlust ./ max(wanderlust);


% start from num_pts points
pts = linspace( min( wanderlust ), max( wanderlust ), num_pts );

% find all points in window_width diameter around each pts
for pt_idx = 1:num_pts
	pt = pts( pt_idx );
	pt_indices{pt_idx} = find( pt - window_width / 2 < wanderlust & wanderlust < pt + window_width / 2 );
end
 
% for each point, calculate mean of each marker around that point
pt_mean_intensity = zeros(num_pts, length(legend_labels));
pt_std_intensity  = zeros(num_pts, length(legend_labels));

if size( data, 2 )>1
    for pt_idx = 1:num_pts
        for marker_idx = 1:size( data, 2 )
            d = data( pt_indices{pt_idx}, marker_idx );
            pt_mean_intensity( pt_idx, marker_idx ) = mean( d );
            pt_std_intensity( pt_idx, marker_idx ) = std( d );
        end
    end
else % calc by grouping
    for pt_idx = 1:num_pts
        d = data( pt_indices{pt_idx});

        pt_mean_intensity( pt_idx, end ) = mean( d );
        pt_std_intensity( pt_idx, end ) = std( d );
        
        % if more than one grouping is specified, calc the mean for each
        if (length(group_ids) > 1)
            groupings = v( pt_indices{pt_idx} );

            pt_mean_intensity( pt_idx, 1:length(group_ids) ) = arrayfun(@(k) mean( d(groupings==k)), group_ids);% ./ mean( d );
            pt_std_intensity( pt_idx, 1:length(group_ids) ) = arrayfun(@(k) std( d(groupings==k)), group_ids);% ./ std( d );
        end
    end
end

% plot
colors = distinguishable_colors( size( pt_mean_intensity, 2 ) );

% shift so the generic blue color is last 
colors = [colors(2:end, :); colors(1, :)];

hold on;
cols = 1:size(pt_mean_intensity,2);
for col = cols
	plot( pts, pt_mean_intensity( :, col ), 'LineWidth', 2, 'Color', colors( col, : ));
end
hold off;


% embellish
% axis( [ 0 1 0 1 ] );
view(2);
legend( legend_labels, 'Location', 'NorthEastOutside', 'Interpreter', 'None' );
graphLabels( '', 'Trajectory', 'Marker intensity (normalized)' );


end
