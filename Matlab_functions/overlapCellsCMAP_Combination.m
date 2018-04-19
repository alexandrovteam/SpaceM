function dummy = overlapCellsCMAP_Combination(alpha, lum_f,cells_p, cmap_p, save_p)
mask = imread(cmap_p);
cell_im = imread(cells_p);

srgb2lab = makecform('srgb2lab');
lab2srgb = makecform('lab2srgb');

shadow_lab = applycform(cell_im, srgb2lab); % convert to L*a*b*

% the values of luminosity can span a range from 0 to 100; scale them
% to [0 1] range (appropriate for MATLAB(R) intensity images of class double) 
% before applying the three contrast enhancement techniques
max_luminosity = lum_f;
L = shadow_lab(:,:,1)/max_luminosity;

% replace the luminosity layer with the processed data and then convert
% the image back to the RGB colorspace
shadow_imadjust = shadow_lab;
shadow_imadjust(:,:,1) = imadjust(L)*max_luminosity;
shadow_imadjust = applycform(shadow_imadjust, lab2srgb);
% 
% shadow_histeq = shadow_lab;
% shadow_histeq(:,:,1) = histeq(L)*max_luminosity;
% shadow_histeq = applycform(shadow_histeq, lab2srgb);
% 
% shadow_adapthisteq = shadow_lab;
% shadow_adapthisteq(:,:,1) = adapthisteq(L)*max_luminosity;
% shadow_adapthisteq = applycform(shadow_adapthisteq, lab2srgb);

% figure, imshow(cell_im);
% title('Original');
% figure, imshow(shadow_imadjust);
% title('Imadjust');
% figure, imshow(shadow_histeq);
% title('Histeq');
% figure, imshow(shadow_adapthisteq);
% title('Adapthisteq');


% alpha = 0.8;
C = alpha * mask + (1 - alpha) * shadow_imadjust;
% imshow(C, [])


% C = imfuse(mask,shadow_imadjust,'blend','Scaling','independent');
% imshow(C, [])
imwrite(C, save_p)
dummy=0;
end
