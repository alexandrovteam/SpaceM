function dummy = imAdjQuantiles(quant, im_p, adj_p)
% quant = 0.02;
% im_p = 'D:\Experiments\Hepa_June\LPS_1.3\Analysis\CellProfilerAnalysis\img_t1_z1_c2.tif';
% adj_p = 'D:\Experiments\Hepa_June\LPS_1.3\Analysis\CellProfilerAnalysis\img_t1_z1_c2_adjusted.tif';

img = imread(im_p);
img = double(img);
img_norm = (img - min(img(:)))/((max(img(:)) - min(img(:))));

low_in = quantile(img_norm(:), quant);
high_in = quantile(img_norm(:), 1-quant);

adjusted = imadjust(img_norm, [low_in high_in], []);
% imshow(adjusted, []);
imwrite(adjusted, adj_p);
dummy = 0;
end