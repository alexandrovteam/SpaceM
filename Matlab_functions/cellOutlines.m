function dummy = cellOutlines(FluoBrightfield_p, fluo_window, label_p, save_p)

% labelI = imread( label_p);%'D:\Experiments\Hepa_May\Analysis\il17a_slide4\Analysis\CellProfilerAnalysis\Labelled_cells.tif');
% fluoI = imread(FluoBrightfield_p);%'D:\Experiments\Hepa_May\Analysis\il17a_slide4\Analysis\ili\FLUO_crop_bin1x1.png');
% values = unique(labelI);
% perimAll = zeros(size(labelI));
% textprogressbar('Drawing cell outlines:     ');

for i = 2:numel(values)
    perim = bwperim(labelI == values(i));
    perimAll = perimAll + perim;
    textprogressbar(i*100/numel(values));
end
PAC = perimAll(fluo_window+1:end-fluo_window, ...
    fluo_window+1:end-fluo_window, :); %PerimAllCut - PAC
fluoI_cut = fluoI(fluo_window+1:end-fluo_window, ...
    fluo_window+1:end-fluo_window, :);

CC = im2double(fluoI_cut).*(~imdilate(PAC, strel('disk', 1))); %CellContours - CC
% imshow(CC, [])

imwrite(CC,save_p);%'D:\Experiments\Hepa_May\Analysis\il17a_slide4\Analysis\CellProfilerAnalysis\Contour_cells.png');
textprogressbar(' \nFinished');
dummy = 0;
end
% folder = 'D:\Experiments\20170410_Hela_STM_mCherry_reanalysis\Input\Microscopy\preMALDI\';
% files = dir(folder);
% for i = 1:numel(files)
%     if strcmp(files(i).name(end-2:end), 'tif') == 1
%         red = imread(strcat(folder,files(i).name),1);
%         dia = imread(strcat(folder,files(i).name),2);
%         new = zeros(1608,1608,2);
%         new(:,:,1) = dia;
%         new(:,:,2) = red;
%         t = Tiff(strcat(folder,'inverted\', files(i).name), 'w');
%         t.write(new);
%     end
% end