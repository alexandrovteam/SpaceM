function dummy = marksMaskCheck(mark_img_p, xye_clean2_p, window, marksMask_p, save_p)

% mark_img_p = 'D:\Experiments\Hepa_June\FA_1.1\Analysis\gridFit\marks_check\PHASE_crop_bin1x1_window100.png';
% xye_clean2_p = 'D:\Experiments\Hepa_June\FA_1.1\Analysis\gridFit\xye_clean2.npy';
% window = 100;
% marksMask_p = 'D:\Experiments\Hepa_June\FA_1.1\Analysis\gridFit\marksMask.mat';
% save_p = 'D:\Experiments\Hepa_June\FA_1.1\Analysis\gridFit\marksMask_check.png';
window = double(window);
load(marksMask_p);
img = imread(mark_img_p);
data = readNPY(xye_clean2_p);
% y = data(1, :) - min(data(1, :)) + window;
% x = data(2, :) - min(data(2, :)) + window;
% waitbar(0)

img_r = img;
img_g = img;
img_b = img;

for i = 1:numel(data(1, :))
    indX = out(i).BWx - min(data(1, :)) + window;
    indY = out(i).BWy - min(data(2, :)) + window;
    
    img_r(sub2ind(size(img), indX, indY)) = 1;
%     img_g(sub2ind(size(img), indX, indY)) = 1;
    img_b(sub2ind(size(img), indX, indY)) = 1;
  
%     waitbar(i/numel(x))
end
rgbImage = cat(3, img_r, img_g, img_b);
% imshow(rgbImage, []);
imwrite(rgbImage, save_p);
dummy = 0;