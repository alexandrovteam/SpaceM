function [newMap ] = interpolate_colormap(colMap, num_values)


    n = num_values;
    m = size(colMap,1);
    t0 = linspace(0,1,m)';
    t = linspace(0,1,n)';
    r = interp1(t0,colMap(:,1),t);
    g = interp1(t0,colMap(:,2),t);
    b = interp1(t0,colMap(:,3),t);
    newMap = [r,g,b];


end

