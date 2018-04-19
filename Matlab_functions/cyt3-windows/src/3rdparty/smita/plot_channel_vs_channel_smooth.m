function [smooth_index, smooth_data] =  plot_channel_vs_channel_smooth(data, channelX, channelY, window_size, colorname)
            
            
            
         
            [num_events, n] = size(data);
            smooth_data = zeros(1,num_events-window_size);
            smooth_index = zeros(1,num_events-window_size);
            
            %need to sort
            
            sorted_data = sortrows(data,channelX);
            
            for i=1:num_events-window_size
                smooth_data(i) = median(sorted_data(i:i+window_size,channelY));
                smooth_index(i) = median(sorted_data(i:i+window_size,channelX));
            end
            
            
            plot(smooth_index,smooth_data, colorname);
        
            
            
end
