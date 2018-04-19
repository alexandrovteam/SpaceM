img = imread('D:\Experiments\20170316_LR_DMAN10_att36_10F1_50x50_50x50\Analysis\StitchedMicroscopy\postMALDI_FLR\img_t1_z1_c1');
Xtot = [];
Ytot = [];
imshow(img, []);hold on;

[x, y] = getpts;
Xtot = [Xtot; x];
Ytot = [Ytot; y];
hold on;scatter(Xtot, Ytot, 20, 'g', 'fill');
[x, y] = getpts;

