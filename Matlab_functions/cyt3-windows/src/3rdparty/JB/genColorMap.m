function [map] = genColorMap(colortag, numLevels)
%% generate colormap
%% default, blue (high) vc yellow (low) map, with 10 levels
%% maxcolor and mincolor should be vectors (r, g, b), real numbers [0 1]
% colortag: matrix containing [maxcolor; mincolor] or [max; mid; min]
% numLevels: number of levels

    if nargin < 2
        numLevels = 10;
    end
    if nargin < 1
        colortag = 'by';        
    end
    
    if isstr(colortag)
        [colortag] = findmap(colortag);
        if isscalar(colortag)
            fprintf('cannot find map\n');            
        end
    end
    
    [nc] = size(colortag,1);
    if nc == 3 && mod(numLevels,2) == 0        
        numLevels = numLevels + 1;
    end
    
    map = zeros(numLevels, 3);
    %if there are three colors, compute colormap twice using numsteps
    numsteps = floor((numLevels-1)/(nc-1));
    for i = 1:nc-1
        steps = (colortag(i,:) - colortag(i+1,:)) ./ numsteps;
        
        map(1+(nc-i)*numsteps,:) = colortag(i,:); %maxcolor
        map(1+(nc-i-1)*numsteps,:) = colortag(i+1,:); %mincolor
        for j = 2:numsteps
            map(j+(nc-i-1)*numsteps,:) = map(j-1+(nc-i-1)*numsteps,:) + steps;
        end
    end
end

function [colortag] = findmap(mapname)    
    ncode = length(mapname);
    colortag = zeros(ncode, 3);
    for i = 1:ncode
        switch mapname(i)
            case 'r'
                colortag(i,:) = [1 0 0];
            case 'g'
                colortag(i,:) = [0 1 0];
            case 'b'
                colortag(i,:) = [0 0 1];
            case 'y'
                colortag(i,:) = [1 1 0];
            case 'w'
                colortag(i,:) = [1 1 1];
            case 'k'
                colortag(i,:) = [0 0 0];
            case 'o'
                colortag(i,:) = [251 170 39]/255;
            case 'p'
                colortag(i,:) = [102 0 145]/255;
            otherwise
                colortag = -1;
                return
        end
    end    
                
%     switch mapname
%         case 'by'
%             colortag(1,:) = [0 0 1];
%             colortag(2,:) = [1 1 0];
%         case 'rw'
%             colortag(1,:) = [1 0 0];
%             colortag(2,:) = [1 1 1];
%         case 'rwb'
%             colortag(1,:) = [1 0 0];
%             colortag(2,:) = [1 1 1];
%             colortag(3,:) = [0 0 1];
%         case 'bwr'
%             colortag(1,:) = [0 0 1];
%             colortag(2,:) = [1 1 1];
%             colortag(3,:) = [1 0 0];
%         otherwise
%             colortag = -1;
%     end           
end


        