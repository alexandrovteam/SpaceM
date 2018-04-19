function nanimagesc(X, cmap, varargin)
    %imagesc for X containing NaN
    %X: data matrix, #condition x #observation
    %cmap: colormap    
    %clim: pass to imagesc
    %'rank': rankdata, true or false
    %'zeroc': fix 0 to mid-level color (eg. black), true or false
    %'nakeep','tiekeep','ascend': {true,false}, only for 'rank'=true
    %example nanimagesc(data, genColorMap('rwb',128), 'clim', [-1 1]))
    
    %default
    zerocenter = false;
    ranked = false;
    nakeep = true;
    tiekeep = true;
    ifascend = 'ascend';
    passpara = cell(0,2);
    nancolor = [0.7 0.7 0.7];
    
    i = 1;
    while i < length(varargin)
        switch lower(varargin{i})
            case 'rank'
                ranked = varargin{i+1};
            case 'zeroc'
                zerocenter = varargin{i+1};            
            case 'ascend'
                if varargin{i+1}
                    ifascend = 'ascend';
                else
                    ifascend = 'descend';
                end
            case 'nancolor'
                nancolor = varargin{i+1};
            case {'xlabels','ylabels','clim','nakeep','tiekeep', 'xbins', 'ybins'}
                eval(sprintf('%s = varargin{i+1};',varargin{i}));            
            otherwise
                %error('Unknown option %s.\n',varargin{i});
                passpara(end+1,1:2) = varargin(i:i+1);
        end
        i = i + 2;
    end
    
    if ranked
        X = rankData(X',ifascend,tiekeep,nakeep)';
    end
    
    nclevel = size(cmap,1); %nclevel is an odd number
    cmap(2:end+1,:) = cmap;    
    cmap(1,:) = cmap(2,:); %nancolor; %set NaN to gray, NaN is consider min in imagesc
    
    if exist('clim','var')
        X(X<clim(1)) = clim(1);
        X(X>clim(2)) = clim(2);
        m = clim(1);
        M = clim(2);        
    else
        m = nanmin(nanmin(X));
        M = nanmax(nanmax(X));        
    end
    if isnan(m)
        m = 0;
    end
    if isnan(M)
        M = 1;
    end
    clim_min = m-1;
    clim_max = M;
    if zerocenter 
        if ~(M > 0 && m < 0)
            M = max(abs([M m]))+1;
            m = -max(abs([M m]))-2;
        end
        
        midcolor = (nclevel+1)/2 + 1; %index of mid-level color in colormap
        halfcolor = (nclevel-1)/2;
        if M > -m
            i = ceil(-m/M*halfcolor);
            cmap = [cmap(1,:); cmap(midcolor-i:end,:)];
        else
            i = ceil(-M/m*halfcolor);
            cmap = [cmap(1,:); cmap(2:midcolor+i,:)];
        end
        nclevel = size(cmap,1) - 1;

        step_size = (M-m)/(nclevel-1);
        colorbar_min = m-step_size; % add an extra step for NaN color
        [~, white_idx] = ismember([1 1 1], cmap, 'rows');

        clim_min = -(white_idx-0.5)*step_size;
        clim_max = (size(cmap,1) - (white_idx-0.5))*step_size;
    else
        step_size = (M-m)/(nclevel-1);
        clim_min = m-step_size;
        clim_max = M;
    end
    
	if exist('xbins','var') 
       
        imagesc(xbins, ybins, X, [clim_min clim_max]);
    else 
        imagesc(X, [clim_min clim_max]);
    end
    
    colormap(cmap);
    if exist('ylabels','var')
        set(gca, 'YTick', 1:length(ylabels), 'YTickLabel',ylabels, 'FontSize',16);
    end
    if exist('xlabels','var')
        set(gca, 'XTick', 1:length(xlabels), 'XTickLabel',xlabels, 'FontSize',18);
    end
    if ~isempty(passpara)
        set(gca, passpara(:,1), passpara(:,2)');
    end
end