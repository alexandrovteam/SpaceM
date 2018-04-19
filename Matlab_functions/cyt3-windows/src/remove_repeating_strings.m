function result = remove_repeating_strings(strings)
    if isempty(strings)
        return;
    end
    
    if numel(strings) == 1
        result = strings;
        return;
    end
    
    str = strings{1};
    for i=1:numel(str)
        TF = strncmpi(str(1:i),strings,i);
        if ~isempty(find(TF==0))
            break;
        end
    end
    
    if i>1
        result = cellfun(@(str)(str(i:end)), strings, 'UniformOutput', false);
    else
        result = strings;
    end
end