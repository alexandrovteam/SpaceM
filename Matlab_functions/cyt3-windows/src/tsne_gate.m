function [ gated_indices ] = tsne_gate(plothandle, mapped_data, ~, type)

    if strcmp(type, 'poly')
        h = impoly(gca);
    elseif strcmp(type, 'ellipse')
        h = imellipse(gca);
    else
        h = imrect(gca);

        disp('Adjust your selection. Double click on node when finished');
        rect = wait(h);
        disp('Processing selection...');

        left = rect(1);
        bottom = rect(2);
        width = rect(3);
        height = rect(4);

        gated_indices = find((mapped_data(:,1)>left) & (mapped_data(:,1)<left+width) & (mapped_data(:,2)>bottom) & (mapped_data(:,2)<bottom+height));
        return;
    end
%     cyt('setStatus', 'Adjust your selection. Double click on node when finished');
    disp('Adjust your selection. Double click on node when finished');
    vert = wait(h);
    disp('Processing selection...');

    gated_indices = find(inpoly(mapped_data, vert));
end