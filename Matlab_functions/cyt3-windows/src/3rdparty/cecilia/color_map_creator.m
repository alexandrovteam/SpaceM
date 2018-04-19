function color_map = color_map_creator(color_scheme)

    %creates color maps, ie N*3 matrices coresponding to RGB colors
    
    color_map = zeros(64,3);
    
    if strcmp(color_scheme, 'yb'),
        %lightyellow-yellow-black-blue-lightblue
        color_map(33:49,:) = [linspace(0,1,17)',linspace(0,1,17)',linspace(0,0,17)']; %black to yellow
        color_map(49:end,:) = [linspace(1,1,16)',linspace(1,1,16)',linspace(0,0.65,16)'];   %yellow to light yellow
        color_map(1:16,:) = [linspace(0.65,0,16)',linspace(0.65,0,16)',linspace(1,1,16)'];
        color_map(16:32,:) = [linspace(0,0,17)',linspace(0,0,17)',linspace(1,0,17)'];

    elseif strcmp(color_scheme, 'rg'),
        %lightred-red-black.green.lightgree colormap
        color_map(33:end,1) = linspace(0,1,32);    %red
        color_map(49:end,2:3) = [linspace(0,0.65,16)',linspace(0,0.65,16)'];    %red
        color_map(1:32,2) = linspace(1,0,32);  %green
        color_map(1:16,[1,3]) = [linspace(0.65,0,16)',linspace(0.65,0,16)'];    %green
        
    elseif strcmp(color_scheme, 'grey_scale_bw'),
        %grey scale color map
        color_map = [linspace(0,1,64)',linspace(0,1,64)',linspace(0,1,64)'];
        
    elseif strcmp(color_scheme, 'grey_scale_wb'),
        %grey scale color map
        color_map = [linspace(1,0,64)',linspace(1,0,64)',linspace(1,0,64)'];
        
    elseif strcmp(color_scheme, 'bw'),
        %black/white colormap, contains only those two colors
        colormap = [1,1,1;0,0,0];

    end
    
end