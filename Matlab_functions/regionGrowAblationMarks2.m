function dummy = regionGrowAblationMarks2(mark_check_p, xye_clean2_p, window, save_p)

mark_check_p = 'D:\Experiments\20171106_Hepa_Nov_DHB_10conditions\FT1\Analysis\gridFit\marks_check\PHASE_crop_bin1x1_window100.png';
% mark_check_p = 'D:\Experiments\20171106_Hepa_Nov_DHB_10conditions\FT1\Analysis\gridFit\marks_check\FLUO_crop_bin1x1_window100.png';
xye_clean2_p = 'D:\Experiments\20171106_Hepa_Nov_DHB_10conditions\FT1\Analysis\gridFit\xye_clean2.npy';
window = 100;
save_p = 'D:\Experiments\20171106_Hepa_Nov_DHB_10conditions\FT1\Analysis\gridFit\marksMask.mat';

cIM = double(imread(mark_check_p));
cIM = (cIM - min(cIM(:))) / (max(cIM(:)) - min(cIM(:)));
cIM = cIM.*-1;
data = readNPY(xye_clean2_p);
y = data(1, :) - min(data(1, :)) + window;
x = data(2, :) - min(data(2, :)) + window;
initPos = round([y;x])';

% figure()
% imshow(cIM, []); hold on;
% 
% ff = fft2(cIM); % Take Fourier Transform 2D
% F1 = 20*log(abs(fftshift(ff)));
% figure();
% imagesc(F1)
% colormap(bone)
% 
% mask1 = (F1- imgaussfilt(F1, 30)); % Blur FT image and subtract it from the original: Highlights the peaks corresponding to high repetition
% mask1(mask1 < 0) = 0; %zoroes values inferior to zero
% imshow(mask1, [])
% mask2 = imdilate(mask1>mean(mask1(:)) + 3.5*std(mask1(:)), strel('disk',2));
% imshow(mask2)
% ff_masked = fftshift(mask2).*ff;
% F2 = 20*log(abs(fftshift(ff_masked)));
% figure();
% imagesc(F2)
% colormap(bone)
% 
% freq_up = 0.75;
% freq_down = 0;
% [N,M] = size(cIM); %[height, width]
% %Sampling intervals
% dx = 1;
% % dy = 1;
% %Characteristic wavelengths
% KX0 = (mod(1/2 + (0:(M-1))/M, 1) - 1/2);
% KX1 = KX0 * (2*pi/dx);
% KY0 = (mod(1/2 + (0:(N-1))/N, 1) - 1/2);
% KY1 = KY0 * (2*pi/dx);
% [KX,KY] = meshgrid(KX1,KY1);
% %Filter formulation
% lpf = (KX.*KX + KY.*KY < freq_up^2);
% lpf2 = (KX.*KX + KY.*KY < freq_down^2);
% mix = lpf.*~lpf2;
% figure();imshow(fftshift(mix).*F2,[])
% colormap(bone)
% 
% rec = real(ifft2(mix.*ff_masked));
% figure(); imshow(rec, [])
% rec = (rec - min(rec(:))) ./ (max(rec(:)) - min(rec(:)));
% % rec = (rec)./max(rec(:));
% rec_adj = imadjust(rec);

rec_adj = cIM;

figure();
imshow(rec_adj);hold on;

range = getrangefromclass(rec_adj);
% waitbar(0);

for i = 1:numel(x)
    out(i).BWx = [];
    out(i).BWy = [];
end

% print('Par Loop start')
% warning('off','all')
% parfor i = 1:numel(x)

for i = 1530:15350
    cw = 30; %cut window
    imCut = rec_adj(initPos(i,1)-cw:initPos(i,1)+cw, initPos(i,2)-cw:initPos(i,2)+cw);
%     imshow(imCut, []);
%     scatter(initPos(i,2),initPos(i,1), 40, 'g', 'fill'); hold on;
% end

    thresh = graythresh(imCut);
    thresh_bw = thresh*range(2);
%     thresh_bw = thresh_bw - thresh_bw*0.5;
    
    thresh_bw = mean(imCut(:)) + std(double(imCut(:)));
    if thresh_bw < 0
        thresh_bw = 0;
    end
    
    indS = find(imCut>=quantile(imCut(:), 0.99));
    if numel(indS) > 0
%         
%         pick = round(numel(indS)/2);
%         [ind1, ind2] = ind2sub([cw*2+1,cw*2+1], indS(pick));
        
        [indSx, indSy] = ind2sub([cw*2+1,cw*2+1], indS);
        indKNN = knnsearch([indSx, indSy], [cw cw]);
        RGpos = [indSx(indKNN) + initPos(i,1)-cw, indSy(indKNN) + initPos(i,2)-cw];
        
%         imshow(imCut>=quantile(imCut(:), 0.99))
%         imshow(imCut,[])
%         hold on;scatter(30,30, 120, 'r', 'fill');
%         scatter(30, 33, 120, 'g', 'fill');
%         
        
%         imshow(rec_adj, []);
%         hold on;
%         scatter(initPos(i,2), initPos(i,1), 10, 'fill', 'r');
%         scatter(ind2 + initPos(i,2)-cw, ind1 + initPos(i,1)-cw, 10, 'fill', 'g');
        
%         rec_adj(RGpos(1), RGpos(2));
        
%         RGpos = [initPos(i,1), initPos(i,2)];
        
        [poly, BW] = regionGrowing(rec_adj, RGpos, thresh_bw, 15);
        [BWx,BWy] = ind2sub(size(rec_adj),find(BW>0));
        scatter(BWy, BWx, 20, 'b', 'fill'); hold on;
        scatter(initPos(i,2),initPos(i,1), 40, 'r', 'fill'); hold on;
        scatter(RGpos(2),RGpos(1), 40, 'g', 'fill'); hold on;

        out(i).BWx = BWx +  min(data(1, :)) - window;
        out(i).BWy = BWy +  min(data(2, :)) - window;
%      
    else strcat('Issue with ablation mark number_', string(i))
    end
    
    % waitbar(i/numel(x))
   

    hold on;
    fprintf('%.0f/%.0f\n', i, numel(x));
%     print(strcat(num2str(i), '/', num2str(numel(x))));     
end
% warning('on','all')
save(save_p, 'out')  
dummy = 0;
end


