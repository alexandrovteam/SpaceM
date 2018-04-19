function centroids = spotFinder(path, matrix_name, fig)
% path = 'D:\Experiments\Analysis_test5\Analysis\StitchedMicroscopy\postMALDI_FLR\Seq0000_XY090.tif';
path = 'Z:\rappez\20170622_Hepa_June_DAN_Untreated_FA_LPS_TNFa\Data\transformation_1\LPS_3.1_SELECTED\Analysis\StitchedMicroscopy\postMALDI_FLR\img_XY066.tif'
% Initial guess of params:
%     freq = 0.20, Rmin = 8, Rmax = 20
% freq = 0.192;
img = imread(path,1);
if numel(size(img)) > 2
    img = img(:,:,1);
end

% figure, imshow(imbinarize(img))
ff = fft2(img); % Take Fourier Transform 2D
F = 20*log(abs(fftshift(ff)));
% figure();
% imagesc(F)
% colormap(bone)

%NEW APPROACH: 
mask = (F- imgaussfilt(F, 15)); % Blur FT image and subtract it from the original: Highlights the peaks corresponding to high repetition
% 
mask(mask < 0) = 0; %zoroes values inferior to zero
rec = real(ifft2(fftshift(imdilate(mask>mean(mask(:)) + 6*std(mask(:)), strel('disk',5))).*ff)); %threshold it, dilate and mask the FT
imshow(rec, [])
% % 
%


        % OLD APPROACH: CREATE DOGHNUT SHAPE TO FILTER FREQUENCIES
        % close all
        freq_up = 0.7;
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
        figure();imagesc(20*log(abs(fftshift(mix.*ff))));colormap(bone);
%         Filter Application
        rec = ifft2(mix.*ff);
        
        

rec = (rec)./max(rec(:));
figure()
imshow(rec, [])

if strcmp(matrix_name,'DMAN') == 1
    
BW = rec > 0.5;
% img = img./max(img(:));
% close all
% se = strel('disk',25);
% tophatFiltered = imtophat(imadjust(rec),se);
%     figure,
end
if strcmp(matrix_name, 'DHB') == 1
rec_adj = imadjust(rec);
% imshow(rec_adj)
% rec_topH = imtophat(rec_adj,strel('disk', 8));
% imshow(rec_topH, [])
BW = imbinarize(rec_adj);
end
% imshow(BW);
BW_c = imclose(BW, strel('disk', 4));
% imshow(BW_c);
BW_o = imopen(BW_c, strel('disk', 7));
% imshow(BW_o);

if bwarea(BW_o) < 1e6
    s = regionprops(BW_o,'centroid');
    centroids = cat(1, s.Centroid);
    
%     a = regionprops(BW_o,'Area');
%     areas = cat(1, a.Area);
%     ind = find(areas> mean(areas+2*std(areas)));
%     ind = find(areas > 400);
%     centroids(ind, :) = [];

else
    centroids = [];
end
% numel(centroids)
if fig 
    figure;
    imshow(img, []); hold on;
    if ~isempty(centroids)
        scatter(centroids(:,1),centroids(:,2), 70, 'y', 'fill');
%     scatter(centroids(ind,1),centroids(ind,2), 70, 'r', 'fill');
    end
    figure;
        subplot(1,2,1);
    imshow(BW_o, []); hold on;
    if ~isempty(centroids)
    scatter(centroids(:,1),centroids(:,2), 70, 'g', 'fill');
%     scatter(centroids(ind,1),centroids(ind,2), 70, 'r', 'fill');
    end
        subplot(1,2,2);
    imshow(rec, []); hold on;
    if ~isempty(centroids)
    scatter(centroids(:,1),centroids(:,2), 70, 'g', 'fill');
%     scatter(centroids(ind,1),centroids(ind,2), 70, 'r', 'fill');
    end
end
% ind = find(areas> mean(areas+2*std(areas)));

%ATTEMPT to watershed --> see watershed documentation
%     D = bwdist(~BW_o);
%     D = -D;
%     D(~BW_o) = -Inf;
%     L = watershed(D);
%     imshow(BW_o, [])
%     imshow(D, [])
%     imshow(L, [])



%      imlabel(BW)
% SEPARATE OBJECTS FROM BW
%     s = regionprops(BW,'centroid');
%     centroids = cat(1, s.Centroid);

%     % FIND CIRCLES ON BW
%     close all
%     [centers, radii, metric] = imfindcircles(BW,[round(Rmin) round(Rmax)]);
end