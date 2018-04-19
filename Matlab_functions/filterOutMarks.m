function datas_2 = filterOutMarks(npy_path, img_path)
% datas = readNPY('D:\Experiments\20171106_Hepa_Nov\FL4\Analysis\gridFit\ablation_marks_XY_2ndFT.npy');
% img = imread('D:\Experiments\20171106_Hepa_Nov\FL4\Analysis\gridFit\blue_window200.png');
datas = readNPY(npy_path);
img = imread(img_path);

f = figure();

imshow(img, []); hold on;
scatter(datas(:,1),datas(:,2), 70, 'g');
datas_2 = datas;

while ishandle(f)
    try
        akZoom(gca)
        [pind,xs,ys] = selectdata('sel','br','action','delete', 'BrushSize', 0.015) ;
        inds = round(pind{1,1});
        datas_2(inds,:) = [];
        
    catch waitforbuttonpress
    end
end
scatter(datas_2(:,1), datas_2(:,2))
end
