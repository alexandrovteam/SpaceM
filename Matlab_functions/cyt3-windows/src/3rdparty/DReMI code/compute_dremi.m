 function [dremi, pvalue, samples_x, samples_y] = compute_dremi(data, noise_threshold , varargin)
            
            
            
            
            
             compute_pvalue = 0;
             set_maxy  = 0;
             num_permutations = 0; 
             max_yval = 0; 
             num_slices = 8; 
             
             
             for i=1:length(varargin)
                 
                if(strcmp(varargin{i}, 'compute_pvalue'))
                    compute_pvalue = 1;
                    num_permutations = varargin{i+1};
                end
                if(strcmp(varargin{i}, 'maxy'))
                    set_maxy = 1;
                    max_yval = varargin{i+1};
                end
                
                if(strcmp(varargin{i}, 'num_slices'))
                    num_slices = varargin{i+1};
                end
                 
                
             end

             [points_x, points_y, point_weights, ~, normalized_density, xaxis,yaxis] = pairwise_visualize(data,'non_averaged_pts', noise_threshold,'no_plot');
            
                
              total_slice_samples = sum(point_weights) * 1000;
              samples_x=[];
              samples_y=[];
             
             for i = 1:length(points_x)
                 
                 num_repeats = point_weights(i) * 1000;
                 new_samples_x = ones(length(num_repeats),1).*points_x(i);
                 new_samples_y = ones(length(num_repeats),1).*points_y(i);
        
                 samples_x = [samples_x; new_samples_x];
                 samples_y = [samples_y; new_samples_y];
             end
                
             
             
            data = [samples_x samples_y];
            size_data = size(data)
            %[dremi, entropy_y, cond_entropy_y] = dremi_resampled(data, 8, 8);
            minx = min(data(:,1));
            miny = min(data(:,2));
            maxx = max(data(:,1));
            maxy = max(data(:,2));
            
            if(set_maxy==1)
                maxy = maxy_val;
            end
            
            dremi = delta_entropyreweight_rawdata(data, minx, miny, maxx, maxy, num_slices,num_slices);
            pvalue = 0;
            if(num_permutations == 0)
                compute_pvalue = 0; 
            end
            
            if(compute_pvalue ==1)
                
                dremi_random = zeros(1,num_permutations);
                total_perms_performed = num_permutations;
                num_greater = 0; 
                for i=1:num_permutations
                     dremi_random(i) = delta_entropyreweight_rawdata(data, minx, miny, maxx, maxy, num_slices,num_slices, 'permute_y');
                    
        
                    if(dremi_random(i)>dremi)
                        num_greater = num_greater+1;
                        if(num_greater==10)
                            total_perms_performed = i;
                            break;
                        end
                    end
        
                end
  
                above_dremi = find(dremi_random(1:total_perms_performed)>dremi);
                if(total_perms_performed<num_permutations)
                    pvalue = 10/total_perms_performed;
                else
                    pvalue = (length(above_dremi)+1)/(num_permutations+1);
                end
            end
            
            
   
            
            
        end