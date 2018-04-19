function calc_corr_wieghted_average_on_rob_anal


[filename, pathname, ~] = uigetfile('*.mat', 'Load robustness analysis saved info');
if isequal(filename,0) || isequal(pathname,0)
    return;
end

load([pathname filename]);
if (isempty(gate_source_corr) || isempty(gate_source_count))
    uiwait(msgbox('Load failed.','Check matlab standard output for error','modal'));    
end

if size(gate_source_corr, 1) == size(gate_source_corr, 2)
    gate_source_corr = gate_source_corr(:, 1)';
    gate_source_count = gate_source_count(:, 1)';
end


nanIdx = isnan(gate_source_corr);

% -- compute and plot pdist corrolation on the 2D tsne channels per source community
sample_shared = gates{end, 2};
gate_sources = unique(sessionData(sample_shared, end-4))';
gate_sources = setdiff(gate_sources, find(nanIdx));

h = figure;
hold on;
colors = distinguishable_colors(max(gate_sources));
gate_names = remove_repeating_strings(gates(gate_sources, 1));
for gate_source_idx = gate_sources
    plot( gate_source_count(gate_source_idx), gate_source_corr(gate_source_idx), '.','MarkerSize', 26, 'Color', colors(gate_source_idx,:) );
end
legend(gate_names, 'Location', 'NorthEastOutside');

% -- magnify legend marker 
l = findobj(gcf,'tag','legend');
a=get(l,'children');

% a(i%3) corresponds to the marker object
for k=1:3:size(a, 1)
    set(a(k),'markersize',24); 
end


% specify plot parameters
axis([0 1000 0 1]);
box on;
xlabel( 'size of manually gated cell type' );
ylabel( 'correlation (pearson)' );
% -- weighted average computation
gate_source_count(nanIdx) = [];
gate_source_corr(nanIdx) = [];

% calc weighted mean and place it on the title
weighted_mean = sum(gate_source_count.*gate_source_corr)./sum(gate_source_count)
title(sprintf('Weighted mean: %g',weighted_mean));
print(h, '-dpng', [pathname 'plot_' filename]);
hold off;

end