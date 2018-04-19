function plot_as_function(x, Y, varargin)
% plots Y (NxM) as a function of vector x (N length).  
% an average of points in Y is computed for 'num_locs' along vector x.
%
% 'labels': a string cell array to enter for the legend
%
% specify 'avg_type': default 'gaussian_var'
%
% 'sliding': a fixed averaged is computed within a sliding window. takes
% into account points within two windows to each side. so let's say 100
% locations would mean each location is averaging points from two windows to
% the right and to the left
%
% 'linear': same as sliding the each point's contibution dimishes linearly
% 
% 'squared': same as sliding the each point's contibution dimishes
% polynomialy
%
% 'gaussian': all points contibute to the weighted average but beyond 2
% standard deviations the contibution is insignificant.
%
% 'gaussian_var': same as gaussian but the standard deviation is variable
% as a function of the normalized density of the points along x. so in
% saprser locations on x, we use a larger standard deviation to let farther
% points contribute.
%
%
% 'show_error': true or false, shows error bars along the curve, default
% 'false'
% 
% Michelle Tadmor, Columbia University, 2013-2015

Markers={'-', ':', ':'};

clear persistent_ksdensity;
legend_flag = false;
avg_type = 'gaussian';
num_locs = 100;
show_error = false;
normalize = true;
rank = false;
smoothness_factor = 0.5;
svGolay = true;
control_density = false;
matPatchColors = [0.75 0.75 0.75; 0.6 0.6 0.6; 0 0 1];
branch = zeros(1, numel(x));
branchY = zeros(1, numel(x));
Y_scale = zeros(1, numel(x));
changes_view = false;
merge_similar = true;
highlight_branch = 0;

for i=1:length(varargin)-1
    if(strcmp(varargin{i},'num_locs'))
        num_locs = varargin{i+1};
    elseif(strcmp(varargin{i},'avg_type'))
        avg_type = varargin{i+1};
    elseif(strcmp(varargin{i},'labels'))
        legend_flag = true;
        labels = varargin{i+1};
    elseif(strcmp(varargin{i},'show_error'))
        show_error = varargin{i+1};  
    elseif(strcmp(varargin{i},'changes_view'))
        changes_view = varargin{i+1};  
    elseif(strcmp(varargin{i},'normalize'))
        normalize = varargin{i+1};
    elseif(strcmp(varargin{i},'rank'))
        rank = varargin{i+1};
        if rank
            control_density = false;
        end
    elseif(strcmp(varargin{i},'svGolay'))
        svGolay = varargin{i+1};
    elseif(strcmp(varargin{i},'smooth'))
        smoothness_factor = varargin{i+1};
    elseif(strcmp(varargin{i},'branch'))
        branch = varargin{i+1};
    elseif(strcmp(varargin{i},'branchY'))
        branchY = varargin{i+1};
        Y_scale = branchY-median(branchY);
        Y_scale(Y_scale<0) = Y_scale(Y_scale<0)./prctile(abs((Y_scale(Y_scale<0))), 99.8);
        Y_scale(Y_scale>0) = Y_scale(Y_scale>0)./prctile(Y_scale(Y_scale>0), 99.8);
        Y_scale(Y_scale<-.3) = -1;        
        Y_scale(Y_scale>.3) = 1;
        
    elseif (strcmp(varargin{i},'highlight'))
        highlight_branch = varargin{i+1};
    end
end

if (rank)
    x = tiedrank(x);
end    

real_weights = zeros(num_locs, length(x));
weights_win = zeros(num_locs, length(x));

tic;
% compute a weight for each value (data point), at each plot location
for i=1:num_locs
    real_weights(i, :) = compute_weights(x, (i/num_locs)*range(x)+min(x), avg_type, smoothness_factor);
    weights_win(i, :)  = compute_weights(x, (i/num_locs)*range(x)+min(x), 'sliding', .1);
end
fprintf('weights computed: %gs\n', toc);

tic;

% Compute weighted averages at each location
% real_weights=robustweighing(real_weights);  % TODO
y_vals_all = real_weights*Y./repmat(sum(real_weights, 2), 1, size(Y, 2));
fprintf('Eighted marker values computed: %gs\n', toc);
    
Y_vals_branches = cell(1,2);
Y_vals_raws     = cell(1,2);

X = linspace(min(x), max(x), num_locs);

if any(Y_scale) 
    % Compute both sides of wine glass
    for bri=1:2
        Y_scale = (-1^(bri-1)) * Y_scale;

        % compute new weights ignore points with a negative Branch
        % association score (BAS of another branch)
        weights{bri}=real_weights.*repmat(1-(max(Y_scale(:)', 0)).^0.9, num_locs, 1);
        fprintf('correcting weights for transitioning: %gs\n', toc);

        % gradually ignore trunk points near the end of the trajectory
        for nl=ceil(num_locs/2):num_locs
            correct_bias = ((num_locs-nl)/ceil(num_locs/2));
            weights{bri}(nl, abs(Y_scale)<0.3) = weights{bri}(nl, abs(Y_scale)<0.3)*correct_bias;
            weights_win(nl, abs(Y_scale)<0.3) = weights_win(nl, abs(Y_scale)<0.3)*correct_bias;
        end

        % correction for the 'short branch'
        branch_sparsity = find(sum(weights_win(:, Y_scale<=-0.2), 2)' > 6);
        if (branch_sparsity(end) < num_locs)      
            weights{bri}(branch_sparsity(end):num_locs, :) =...
                repmat(weights{bri}(branch_sparsity(end), :), num_locs-branch_sparsity(end)+1, 1);
        end

        % Compute weighted averages at each location
        Y_vals_raws{bri} = weights{bri}*Y./repmat(sum(weights{bri}, 2), 1, size(Y, 2));
    end
else
    % Compute weighted averages at each location
    Y_vals_raws{1} = real_weights*Y./repmat(sum(real_weights, 2), 1, size(Y, 2));
    
    % pretend like there are two identical branches to simplify the rest of
    % the code
    Y_vals_raws{2} = Y_vals_raws{1};
    weights{1}=real_weights;
	weights{2}=real_weights;
end

tic

% we want to normalize to [0 1] branches as well as trunk
both_branches = [Y_vals_raws{1};Y_vals_raws{2}];

mins = prctile(both_branches, 0, 1);
rngs = prctile(bsxfun(@minus, both_branches, mins), 100, 1); 

if normalize  
    y_vals_all=translate(y_vals_all, mins, rngs);
    Y_vals_branches{1} = translate(Y_vals_raws{1}, mins, rngs);
    Y_vals_branches{2} = translate(Y_vals_raws{2}, mins, rngs);
else
    Y_vals_branches = Y_vals_raws;
end
fprintf('Normalization values computed: %gs\n', toc);

Y_vals_main = Y_vals_branches{1};
Y_vals_brnch = Y_vals_branches{2};

% branch line - we are selective on the plotting
% some markers stay the same on both branches so we want to keep them
% together
if (merge_similar)
    % look at the differences between the two branches
    diffs = abs(Y_vals_brnch - Y_vals_main);
    yrange = max([max(Y_vals_brnch), max(Y_vals_main)]) - min([min(Y_vals_brnch), max(Y_vals_main)]);
    for mi=1:size(Y_vals_brnch, 2)
        span = 10;
        z = smooth(diffs(:, mi), span);
        branch_locations = find(z>.06*yrange);
      
        if numel(branch_locations) < 10
            Y_vals_brnch(:,mi) = nan;
            
        elseif branch_locations(1) > span
            Y_vals_brnch(1:(branch_locations(1)-span), mi) = nan;
            a = branch_locations(1)-span+1;
            b = branch_locations(1);
            ws = (1:span)'./span;
            ws = ws-ws(1);
            Y_vals_brnch(a:b,mi) = (1-ws).*y_vals_all(a:b, mi) + ws.*Y_vals_brnch(a:b, mi);
            Y_vals_main(a:b, mi) = (1-ws).*y_vals_all(a:b, mi) + ws.*Y_vals_main(a:b, mi);
%             Y_vals_main((b-5):(b+5), mi) = smooth(Y_vals_main((b-5):(b+5), mi),3);
        end
    end
end

Y_vals_main(isnan(Y_vals_brnch)) = y_vals_all(isnan(Y_vals_brnch));

if (svGolay)  
    for col=1:size(Y_vals_brnch, 2)
        Y_vals_main(:, col) = smooth(X, Y_vals_main(:, col),sqrt(num_locs*2), 'sgolay');
    end          
end

Y_vals_branches{1} = Y_vals_main;
Y_vals_branches{2} = Y_vals_brnch;

% plot branch ony if a branchY was given
plot_branch=1;
if any(branchY) 
    plot_branch=2; 
end

% iterate for plotting
sz = ones(1,2);
if highlight_branch
    sz(highlight_branch) = 2;
end
    
for bri=1:plot_branch
    matColors = distinguishable_colors(size(Y, 2));
%     matColors = jet(size(Y, 2));
    set(gca, 'ColorOrder', matColors);
    set(0, 'DefaultAxesColorOrder', matColors);

    % change marker selection
    marker = Markers{bri};
    Y_vals = Y_vals_branches{bri};
    
    % show how\when the markers change over time using a heatmap
    if (changes_view)
        
        Y_vals(2:end, :) = diff(Y_vals);
        Y_vals(1, :) = Y_vals(2,:);
        
        axis auto;
        box on;
        hold off;
        xticks = 10;
        frmt = '%1.1f';
        if (max(x) > 99)
            xticks = 5;
            frmt = '%1.0f';
        end
        xlbls = cellfun(@(n)num2str(n, frmt), num2cell(linspace(min(x), max(x), xticks)), 'UniformOutput',false);
        nanimagesc(Y_vals', genColorMap('rwb', 20), 'ylabels', labels, 'zeroc', 1);
        set(gca, 'XTick', linspace(1, numel(X), xticks), 'XTickLabel', xlbls);
        set(gca,'TickLabelInterpreter','none')
        title('Wishbone trajectory' );
        
        colorbar;

        return;
    end
    
    plot(X, Y_vals(:, 1),marker,...
         'LineWidth', 3*sz(bri),...
         'markersize', 4*sz(bri),...
         'Color', matColors(1, :)); 
    hold on;
    if (size(Y, 2)> 1)
        for col=2:size(Y, 2)
            hold on;
            plot(X, Y_vals(:, col),marker,...
             'LineWidth', 3*sz(bri),...
             'markersize', 4*sz(bri),...
             'Color', matColors(col, :));        
        end
    else    
%         hold on;
%         if (normalize)
%             scatter(x, translate(Y, mins, rngs), '.b');
%         else
%             scatter(x, Y, '.b');
%         end
    end
    
    % plot std deviation
    if show_error 
        % compute variace along X (symmetrically)
        for i=1:num_locs       
            
            % symmetrical
            if (normalize)
                Ys = translate(Y, mins, rngs);
            else
                Ys = Y;
            end
            
            Y_errs = bsxfun(@minus,Ys,(Y_vals(i, :)));

            M = sum(weights{bri}(i, :)~=0);
            s = (M-1)/M;
            w_sum = sum(weights{bri}(i, :));

            Y_valerrs(i, :) = sqrt((weights{bri}(i, :)*((Y_errs).^2))/(s*w_sum));
        end

%         if (normalize)
%             Y_valerrs = bsxfun(@rdivide,Y_valerrs,rngs);
%         end

        % plot the variance as a pretty translucent cloud around line
        for yi=1:size(Y, 2)

            % plot light grey background first for variance
            X_fill = [X, fliplr(X)];
            
            Y_i = Y_vals(:, yi)';
            Y_i_err = Y_valerrs(:, yi)';
            
            Y_fill = [Y_i-.5*Y_i_err, fliplr(Y_i+.5*Y_i_err)];
            
            fill( X_fill(~isnan(Y_fill)), Y_fill(~isnan(Y_fill)),...
                matColors(yi,:),'linestyle','none','facealpha',.5);
        end
    end

    % Hold on if plotting more branch lines
	hold on;
end

if normalize && ~changes_view
    ylim([0 1]);
else
    ylim([min(mins) max(mins+rngs)]);
end

    % show density histogram under the plot to show the concentration 
    if false %~(rank || control_density)    
        try
            dens = sum(weights_win, 2)';
            ca = axis;
            hold on;
            y_range = ca(4)-ca(3);
            y_buffer = y_range*.04;
        %     ylim([ca(3)-.1 ca(4)]);
            imagesc(X, [ca(3)-y_buffer ca(3)-(.5*y_buffer)], Y_dens{1}, [0, max(Y_dens{1})]);
            hold on;
%             imagesc(X, [ca(3)-y_buffer ca(3)-(.5*y_buffer)], Y_dens{2}, [0, max(Y_dens{2})]);   
            colorbar;
            axis([ca(1:2) ca(3)-y_buffer ca(4)]);
        catch e
            disp(getReport(e,'extended'));
        end
    end

if (legend_flag)
    l=legend(labels, 'Interpreter', 'none');
%     set(l, 'Location','NorthEastOutside');
end

% if (check==0)
%     % saving into file
%     try
%     save(cachefilename,'mapMat');
%     catch
%         fprintf('error caching weights in %s', cachefilename);
%     end
% end

end

function weights = compute_weights(points, loc, type, factor)
    
    range = quantile(points, .98) - quantile(points, .02);
    min_std_dev = factor*.1*range; % minimum std_dev for dense regions
    max_std_dev = .19*range; % max std_dev for sparse regions  
    linear_slope = 10/range;
    
    if strcmpi('sliding', type) %set '1's on the indices in the windows 
        weights = (points < (loc + 2*(min_std_dev))) & ...
                  (points > (loc - 2*(min_std_dev)));
    
    elseif strcmpi('linear', type)
        weights = 1 - linear_slope*(abs(points - loc));
        weights(weights<0) = 0;
        weights = weights/sum(weights);       
    
    elseif strcmpi('squared', type)
        weights = 1 - ((linear_slope*(points - loc)).^2);
        weights(weights<0) = 0;
        weights = weights/sum(weights);
             
    elseif strcmpi('gaussian_var', type)
        [f, xi] = persistent_ksdensity(points);
    	d = f(minind(abs(xi-loc)))/max(f);

        std_dev = d*(min_std_dev)+(1-d)*max_std_dev;

        weights = ((2*pi*(std_dev^2))^(-1))*exp(-.5*((points - loc)/std_dev).^2);
    
    else % default is strcmpi('gaussian', type)
        weights = ((2*pi*(min_std_dev)^2)^(-1))*exp(-.5*((points - loc)/min_std_dev).^2);
    end
end

function weights=robustweighing(weights)
    for loc=1:size(weights, 1)
        % no single cell should contribute over 30% of the signal
        to_slash = (weights(loc, :)./sum(weights(loc, :)) > .35);
%         weights(loc, to_slash) = weights(loc, to_slash)./4;
%         
%         weights(loc, :) = weights(loc, :).*(exp(-weights(loc, :)./(max(weights(loc, :)))));
    end
end

function Y=translate(X, min, range)
    Y = bsxfun(@minus, X, min);
    Y = bsxfun(@rdivide, Y, range);

    Y(Y<0) = 0;        
    Y(Y>1) = 1;
end

function [f, xi] = persistent_ksdensity(points)

    persistent f_p;
    persistent xi_p;    

    if isempty(f_p)
        [f_p, xi_p] = ksdensity(points);
    end
    f = f_p;
    xi = xi_p;
end