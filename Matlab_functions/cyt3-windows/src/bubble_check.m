function bubble_check(data, vtime, channel_names, gate_name) 
    %% bubble check
    binsize = 20;    
    nbins = floor(size(data,1)/binsize);
    data = data(1:nbins*binsize, :);
    [~, idx] = sort(data(:, 1));
    try
        h = figure;
        for ch=1:length(channel_names)
                valuesperbin = reshape(data(idx, ch), binsize, nbins);
                means = mean(valuesperbin, 1);
                plot(1:nbins, means , '-');            
                title(channel_names{ch});
                drawnow;
                print(h,'-dpng', validfilename(sprintf('%s %s vs time', gate_name, channel_names{ch})));
        end
    catch e
        close(h);
        disp(getReport(e,'extended'));
        return;
    end
    close(h);
    return;
end