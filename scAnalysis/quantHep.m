function data = quantHep(file_p)

% file_p = 'D:\Experiments\20161219_LR_dhb10_rec2_space50_att32_50x50_hepatocytes\Input\Microscopy\preMALDI\Seq0000_XY140.tif';
% file_p = 'C:\Users\Luca\Desktop\Composite.png';
I1 = imread(file_p, 1);
I2 = imread(file_p, 2);
% I3 = imread(file_p, 3);

IM1 = double(I1);
IM1 = IM1./max(IM1(:));
% imshow(IM1, [])

[centers, radii] = imfindcircles(IM1,[5 50]);
% imshow(IM1, []); hold on;
% scatter(centers(:,1),centers(:,2), 100, 'r', 'fill');

%Reinforce edges on fat label
IM2 = double(I2);
IM2 = IM2./max(IM2(:));
% imshow(imadjust(IM2), []);

% niblack threshold + make the mask more fancy
IM2bw1 = imfill(niblack(IM2, [50 50], 0));
IM2bw2 = imopen(IM2bw1, strel('disk', 3));
IM2bw3 = imerode(imfill(imdilate(IM2bw2, strel('disk', 5))), strel('disk', 6));
bw4 = bwareaopen(IM2bw3, 500); %remove labels with area inferior to 500

bw4_perim = bwperim(bw4);
% overlay1 = imoverlay(IM2, bw4_perim, [.3 1 .3]);
% imshow(overlay1);hold on;
% scatter(centers(:,1), centers(:,2), 100, 'r', 'fill'); %--> markers are outside and inside contours

% find markers inside contours -> red out, green in
in = intersect(find(bw4 == 1), sub2ind(size(bw4), round(centers(:,2)), round(centers(:,1))));
[inx, iny] = ind2sub(size(bw4), in);
X = [inx, iny];
MdlKDT = KDTreeSearcher(X);
[Idx,D] = knnsearch(MdlKDT,X,'K',2);

% overlay3 = imoverlay(imadjust(IM2), bw4_perim, [1 .3 .3]);
% imshow(overlay3);hold on;
% scatter(centers(:,1), centers(:,2), 100, 'r', 'fill');
% scatter(iny,inx, 100, 'g', 'fill');
% C = imfuse(imadjust(IM1), imadjust(IM1) + imadjust(IM2), imadjust(IM1));
% imshow(C, []);

% scatter(iny(find(D(:,2) <= 30)),inx(find(D(:,2) <= 30)), 100, 'g', 'fill' ); 
% hold on; 
% viscircles([iny, inx],ones(288, 1)*15) 

% iny_c = setdiff(iny, iny(find(D(:,2) <= 30),1));
% inx_c = setdiff(inx, inx(find(D(:,2) <= 30),1));

iny(find(D(:,2) <= 30)) = [];
inx(find(D(:,2) <= 30)) = [];


% creaste a mask containing the cell markers -> used later for seeded
% watershed
mask_em = zeros(1608,1608);
mask_em(sub2ind(size(mask_em), inx, iny)) = 1;
me_d = imdilate(mask_em, strel('disk', 14));

meIM2 = double(I2) .* me_d;

s = regionprops(logical(me_d), meIM2, {'Centroid','PixelValues','BoundingBox'});
% imshow(imoverlay(imadjust(IM2), bwperim(me_d), [.3 1 .3]), []);
% title('Mean Intensity of Regions');
% hold on
numObj = numel(s);
data = zeros(numObj, 3);
for k = 1 : numObj
    s(k).MeanIntensity = mean(double(s(k).PixelValues));
%     text(s(k).Centroid(1),s(k).Centroid(2), ...
%         sprintf('%2.1f', s(k).MeanIntensity), ...
%         'EdgeColor','b','Color','g');
    data(k,:) = [s(k).Centroid(1),s(k).Centroid(2),s(k).MeanIntensity];
end
% hold off
% 
% int = [s.MeanIntensity]';

% data = [[s(:,1).Centroid]', [s(:,2).Centroid]', [s.MeanIntensity]'];

% overlay2 = imoverlay(IM2, bw4_perim | mask_em, [.3 1 .3]);
% imshow(overlay2)

% IM2_c = imcomplement(IM1);
% I_mod = imimposemin(IM1, ~bw4 | mask_em);
% L = watershed(I_mod);
% imshow(imoverlay(I_mod,bwperim(L),[1 .3 .3]), []);
% imshow(label2rgb(L))














end