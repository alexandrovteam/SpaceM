function output_txt = clusterLabels(obj,event_obj,index,cluster_size,cellsInCluster)
    % Display information about the choosen cluster
    % obj          Currently not used (empty)
    % event_obj    Handle to event object
    % xydata       Entire data matrix
    % labels       State names identifying matrix row
    % xymean       Ratio of y to x mean (avg. for all obs.)
    % output_txt   Datatip text (string or string cell array)
    % This datacursor callback calculates a deviation from the
    % expected value and displays it, Y, and a label taken
    % from the cell array 'labels'; the data matrix is needed
    % to determine the index of the x-value for looking up the
    % label for that row. X values could be output, but are not.

    pos = get(event_obj,'Position');
    x = pos(1); y = pos(2);
    output_txt = {['Index ' num2str(index)],['Cluster Size: ' num2str(cellsInCluster),' cells (',num2str(cluster_size*100,'%2.2f') '%)']};
end