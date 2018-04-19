function x = treshAndSURF(im1_path, treshold)
% im1_path = 'D:\Experiments\Analysis_test5\Analysis\StitchedMicroscopy\postMALDI_FLR\Seq0000_XY010.tif';
img1 = imread(im1_path, 1);
n_dim = size(size(img1));
if n_dim(2) > 2
    img1 = img1(:,:,1);
end
if isempty(treshold)
    level1 = graythresh(img1);
    BW1 = im2bw(img1,level1);
else
    BW1 = img1 < treshold;
end
SE = strel('disk',20);
% % clearvars('img1');
%imshow(BW1)
IM1 = imclose(BW1,SE);
IM1 = imdilate(IM1,SE);
%imshow(IM1)
% clearvars('BW1');
IM2 = imfill(~IM1,'holes');
% IM2 = imclose(IM2, strel('disk', 100));
a = detectSURFFeatures(IM2);
% clearvars('IM1');
% imshow(img1,[])
% hold on;scatter(a.Location(:,1),a.Location(:,2))
if isempty(a.Location)
    x = [];
else
    x(:,1) = double(a.Location(:,1));
    x(:,2) = double(a.Location(:,2));
end
% print x

% imshow(img1, []);
% hold on;
% scatter(x(:,1), x(:,2), 'g', 'fill');
% img2 = imread(im2_path);
% level2 = graythresh(img2);
% BW2 = im2bw(img2,level2);
% clearvars('img2');
% IM2 = imdilate(BW2,SE);
% clearvars('BW2');
% b = detectSURFFeatures(IM2);
% clearvars('IM2');
end
