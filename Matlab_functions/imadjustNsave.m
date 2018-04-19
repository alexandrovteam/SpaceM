function dummy = imadjustNsave(im_in_path, im_out_path)
im = imread(im_in_path);
im_adj = imadjust(im);
imwrite(im_adj, im_out_path)
dummy = 0;
end