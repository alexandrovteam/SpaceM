function my_plot_dens(vx,vy, norm)
% norm == 0 no cond dens
% norm == 1 cond on vx 
% norm == 2 cond on vy

    if norm==0;
        kde2d_color_hist([vx vy], 'thresh', 0.001);
    elseif norm==1;
        kde2d_color_hist([vx vy], 'thresh', 0.001, 'axis', 2);
    else
        kde2d_color_hist([vx vy], 'thresh', 0.001, 'axis', 1);
    end
end