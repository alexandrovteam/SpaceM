function centroids = spotFinder2(path, layer)
path = 'D:\Experiments\20171106_Hepa_Nov_DHB_10conditions\FT4\Analysis\gridFit\blue_window200.png';
% close('all');
img_i = double(imread(path));
% if numel(size(img_i)) > 2
%     img_i = img_i(:,:,1);
% end

img = img_i(:,:,uint8(layer)).*-1;
img = (img - min(img(:))) ./ (max(img(:)) - min(img(:)));
% figure; imshow(img , [])
low_in = mean(img(:)) + 2*std(img(:)); 
if low_in >=1, low_in = 0.8; end
img_2 = imadjust(img, [low_in, 1]);
figure(); imshow(img_2, [])
ff = fft2(img_2); % Take Fourier Transform 2D
F1 = 20*log(abs(fftshift(ff)));
% figure();
% imagesc(F1)
% colormap(bone)

%NEW APPROACH: 
mask1 = (F1- imgaussfilt(F1, 15)); % Blur FT image and subtract it from the original: Highlights the peaks corresponding to high repetition
mask1(mask1 < 0) = 0; %zoroes values inferior to zero
% imshow(mask1, [])
mask2 = imdilate(mask1>mean(mask1(:)) + 3.5*std(mask1(:)), strel('disk',2));
% imshow(mask2)
ff_masked = fftshift(mask2).*ff;
F2 = 20*log(abs(fftshift(ff_masked)));
% figure();
% imagesc(F2)
% colormap(bone)

% rec = real(ifft2()); %threshold it, dilate and mask the FT
% imshow(rec, [])
% 
% OLD APPROACH: CREATE DOGHNUT SHAPE TO FILTER FREQUENCIES
% close all
freq_up = 0.9;
freq_down = 0;
[N,M] = size(img); %[height, width]
%Sampling intervals
dx = 1;
% dy = 1;
%Characteristic wavelengths
KX0 = (mod(1/2 + (0:(M-1))/M, 1) - 1/2);
KX1 = KX0 * (2*pi/dx);
KY0 = (mod(1/2 + (0:(N-1))/N, 1) - 1/2);
KY1 = KY0 * (2*pi/dx);
[KX,KY] = meshgrid(KX1,KY1);
%Filter formulation
lpf = (KX.*KX + KY.*KY < freq_up^2);
lpf2 = (KX.*KX + KY.*KY < freq_down^2);
mix = lpf.*~lpf2;
% figure();imshow(fftshift(mix).*F2,[])
% colormap(bone)
%         Filter Application
rec = real(ifft2(mix.*ff_masked));
% figure(); imshow(rec, [])

rec = (rec - min(rec(:))) ./ (max(rec(:)) - min(rec(:)));
% rec = (rec)./max(rec(:));
rec_adj = imadjust(rec);
% figure(); imshow(rec_adj, []);

rec_bw = rec_adj > mean(rec_adj(:)) + 2*std(rec_adj(:));
% figure(); imshow(rec_bw);

BW_o = imopen(rec_bw, strel('disk', 3));
% imshow(BW_o);

s = regionprops(BW_o,'centroid');
centroids = cat(1, s.Centroid);
%     
% figure;
% imshow(img_2, []); hold on;
% if ~isempty(centroids)
%     scatter(centroids(:,1),centroids(:,2), 150, 'g', 'LineWidth', 5);
% %     scatter(centroids(ind,1),centroids(ind,2), 70, 'r', 'fill');
% end
% close all
end

