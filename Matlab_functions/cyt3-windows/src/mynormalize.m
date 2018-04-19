function data = mynormalize(data, varargin)
% normalize data according to specified percentile
    percentile = 0;
    for i=1:length(varargin)-1
        if(strcmp(varargin{i},'rows'))
            num_locs = varargin{i+1};
        end
        if (strcmp(varargin{i},'percentile'))
            percentile = varargin{i+1};
        end
    end
    
    if percentile
        fprintf('Normalizing according to the %gth percentile...', percentile);
        data = data-repmat(prctile(data, 100-percentile, 1), size(data,1),1);
        data = data./repmat(prctile((data), percentile, 1),size(data,1),1);
    else
        data_shift = data - repmat(min(data), size(data,1),1);
        data = data_shift./repmat(max(data_shift),size(data,1),1);
    end
    
    data(data > 1) = 1;
    data(data < 0) = 0;
    data(isinf(data)) = 0;
end