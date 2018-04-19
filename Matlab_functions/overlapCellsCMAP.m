function dummy = overlapCellsCMAP(cells_p, cmap_p, save_p)
mask = imread(cmap_p);
cell_im = imread(cells_p);
C = imfuse(mask,cell_im,'blend','Scaling','independent');
imwrite(C, save_p)
dummy=0;
end
