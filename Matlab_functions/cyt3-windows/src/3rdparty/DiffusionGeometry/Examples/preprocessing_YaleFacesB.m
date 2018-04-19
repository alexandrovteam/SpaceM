close all
clear 
clc

% go to the directory where the files are stored
%cd /home/faculty/glchen/downloads/yaleB

I = zeros(480,640,65,9,10);
Labels = zeros(3,65,9,10);
% 480-by-640 images under 65 illumination conditions, in 9 poses, and for 10
% human subjects

for subject = 1:10
    fprintf('\n Subject: %d',subject);
    for pose = 0:8
        
        folder = sprintf('yaleB%.2d_P%.2d',subject,pose);
        eval(['cd ' folder]);
        
        file = [folder '.info'];
        f = fopen(file);
        
        C = textscan(f, '%s');
        C = C{1};

        for illumination = 1:length(C)
            I(:,:, illumination, pose+1, subject) = getpgmraw(C{illumination});
            Labels(:,illumination,pose+1,subject) = [illumination,pose,subject];            
        end
       
        
        ST = fclose(f);
        
        cd ..
        
    end
    
end

%
nImages = 65*9*10;

I = reshape(I, [480, 640, nImages]);
Labels = reshape(Labels,[3,nImages])';

I = reshape(I, [], nImages);

Imean = mean(I,2);

for k = 1:size(I,2);I(:,k)=I(:,k)-Imean;end;

[U,S,V]=pca(I,200);

figure;plot(diag(S))

save YaleB_PCA Labels U S V Imean