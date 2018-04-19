function [centers, radii] = spotFinder(path)
%     path = imread('D:/Experiments/20160810_HeLa_Stm_WT_mCHerry/Sample2/Input/Microscopy/postMALDI/Seq0000_XY077.tif');
    % figure, imshow(imbinarize(img))
    ff = fft2(imbinarize(img)); % Take Fourier Transform 2D
    F = 20*log(abs(fftshift(ff)));
    figure();
    imagesc(F) 
    colormap(bone)

    close all
    K0 = 0.20;
    [N,M] = size(img); %[height, width]
    %Sampling intervals 
    dx = 1; 
    dy = 1; 
    %Characteristic wavelengths 
    KX0 = (mod(1/2 + (0:(M-1))/M, 1) - 1/2); 
    KX1 = KX0 * (2*pi/dx); 
    KY0 = (mod(1/2 + (0:(N-1))/N, 1) - 1/2); 
    KY1 = KY0 * (2*pi/dx); 
    [KX,KY] = meshgrid(KX1,KY1); 
    %Filter formulation 
    lpf = (KX.*KX + KY.*KY < K0^2); 
    %Filter Application 
    rec = ifft2(lpf.*ff);

    rec = rec./max(rec(:));
    % img = img./max(img(:));
    close all
    % se = strel('disk',25);
    % tophatFiltered = imtophat(imadjust(rec),se);
    figure,
    rec_adj = imadjust(rec);
    BW = imbinarize(rec_adj);

    % SEPARATE OBJECTS FROM BW
%     s = regionprops(BW,'centroid');
%     centroids = cat(1, s.Centroid);

    % FIND CIRCLES ON BW
    close all
    [centers, radii, metric] = imfindcircles(BW,[8 20]);
end