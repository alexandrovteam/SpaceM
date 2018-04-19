function [delta, entropy_y, cond_entropy_y, hist_bins, data, edges] = delta_entropyreweight_rawdata(data,minx, miny, maxx, maxy, num_slicesx,  num_slicesy, varargin)
   
    %get the density estimate
    
    sprintf('in delta_entropyreweight');
    
    %channel1 = cdata.name_channel_map(channel1_name);
    %channel2 = cdata.name_channel_map(channel2_name); 
    %data = [cdata.data(:,channel1) cdata.data(:,channel2)];
    %erase the values after the max value
    
%     minx = min(data(:,1));
%     maxx = max(data(:,1));
%     miny = min(data(:,2));
%     maxy = max(data(:,2));
%     
   
%     erase_indices = find(data(:,2)>maxy);
%     data(erase_indices,:) = [];
%     erase_indices = find(data(:,1)>maxx);
%     data(erase_indices,:)=[];
%     erase_indices = find(data(:,2)<miny);
%     data(erase_indices,:) = [];
%     erase_indices = find(data(:,1)<minx);
%     data(erase_indices,:)=[];
%     
  
    fixed_partition = 0;
    
    for i=1:length(varargin)
    
        if(strcmp(varargin{i},'permute_y'))
            %for the purpose of computing a p-value 
            
            p = randperm(length(data));
           
        
            data1 = data(p,1);
            data2 = data(:,2);
            data = [data1 data2];
        end
        
        if(strcmp(varargin{i},'fix_partition'))
            
            edges = varargin{i+1};
            
            fixed_partition = 1;
        end
    end
    
    if(fixed_partition == 0)
        xincrement = (maxx-minx)/num_slicesx;
    	yincrement = (maxy-miny)/num_slicesy;
    
    
       edges_x = [minx];
   
       for i=1:num_slicesx-1;
        
            edges_x = [edges_x edges_x(i)+xincrement];
        
        
       end
    
       edges_x = [edges_x maxx];
    
       edges_y = [miny];

       for i=1:num_slicesy-1;

            edges_y = [edges_y edges_y(i)+yincrement];


       end

       edges_y = [edges_y maxy];

    

        edges = cell(2,1);
        edges{1} = edges_x;
        edges{2} = edges_y;
   
    end
    
    
    hist_bins = transpose(hist3(data,'Edges', edges));
    
    
    
    hist_bin_totals = sum(hist_bins);
    hbin_totals = transpose(hist_bin_totals);
   
    
    %low_cols = find(hist_bin_totals<100);
    %hist_bins(:,low_cols) = [];
    
    %hist_bin_row_totals = sum(hist_bins,2);
    %low_rows = find(hist_bin_row_totals < 20);
    %hist_bins(low_rows,:) = [];
    
    num_rows = size(hist_bins,1);
    num_cols = size(hist_bins,2);
    
    
    
    entropy_y = 0;
   
    
    column_weights = ones(1,num_cols);
    slice_entropies = zeros(1, num_cols);
   
   
    total = sum(column_weights);
    
    
    for i=1:num_rows
        
       p_y = 0; 
       for j=1:num_cols
            column_total = sum(hist_bins(:,j));
           
            
            if(column_total > 0)
                 
                p_y = p_y + (hist_bins(i,j)*(column_weights(j)/column_total));
            end
       end
       
       %at this point I have the appropriately weighted sum of everything
       %in that y slice 
       
       p_y = p_y/total;
        
         if (p_y==0)
                
                p_y = 0;
                log_p_y = 0;
         else
            log_p_y = log2(p_y);
             
         end
        
        entropy_y = entropy_y+(p_y * log_p_y);
    end
   
    entropy_y = -1*entropy_y;
    
    cond_entropy_y  = 0;
    
    for j=1:num_cols
        
       
        slice_entropy_y = 0;
       
        col_sum = sum(hist_bins(:,j));
        
        if(col_sum ==0)
            continue;
        end
        norm_col_hist = hist_bins(:,j)./col_sum;
        
        
        for i=1:num_rows
        
            
            slice_p_y = norm_col_hist(i);
            
            if (slice_p_y==0)
                
                slice_p_y = 0;
                log_slice_p_y = 0;
            else
                slice_p_y;
                log_slice_p_y = log2(slice_p_y);
            end
        
            slice_entropy_y = slice_entropy_y + (slice_p_y * log_slice_p_y);
        end
        
        slice_entropy_y = -1*slice_entropy_y;
        
        slice_entropies(j) = slice_entropy_y;
        
        col_prob = column_weights(j)/sum(column_weights);
        
        cond_entropy_y = cond_entropy_y + (col_prob * slice_entropy_y);
        
    end
    
     slice_entropies;
     
     delta_before_division = entropy_y - cond_entropy_y;
     cond_entropy_y_before_division = cond_entropy_y;
     %normalize this by the size of the grid 
     smaller_dim = min(num_rows, num_cols);
     log_smaller_dim = log2(smaller_dim);
     %this is supposed to be the maximum achievable value on a grid this
     %size
     
     delta = delta_before_division/log_smaller_dim;
     cond_entropy_y = cond_entropy_y/log_smaller_dim;
     
     
    