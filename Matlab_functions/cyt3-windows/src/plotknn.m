function plotknn(D, k)
    IDX = knnsearch(D,D, 'K', k+1);
    for i=1:size(IDX,1)
        for j=2:size(IDX,2)
            pts = [D(i, :);  D(IDX(i, j), :)];
            if (size(D,2) == 2)
                line(pts(:, 1), pts(:,2));
            else
                line(pts(:, 1), pts(:,2), pts(:,3));
            end
        end
    end
end