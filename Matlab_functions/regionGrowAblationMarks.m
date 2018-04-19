function dummy = regionGrowAblationMarks(mark_check_p, xye_clean2_p, window, save_p)

% % mark_check_p = 'D:\Experiments\20171106_Hepa_Nov_DHB_10conditions\FT1\Analysis\gridFit\marks_check\PHASE_crop_bin1x1_window100.png';
% mark_check_p = 'D:\Experiments\20171106_Hepa_Nov_DHB_10conditions\FT1\Analysis\gridFit\marks_check\FLUO_crop_bin1x1_window100.png';
% xye_clean2_p = 'D:\Experiments\20171106_Hepa_Nov_DHB_10conditions\FT1\Analysis\gridFit\xye_clean2.npy';
% window = 100;
% save_p = 'D:\Experiments\20171106_Hepa_Nov_DHB_10conditions\FT1\Analysis\gridFit\marksMask.mat';

cIM0 = double(imread(mark_check_p));
cIM = cIM0.*-1;
cIM = (cIM - min(cIM(:))) / (max(cIM(:)) - min(cIM(:)));
data = readNPY(xye_clean2_p);
y = data(1, :) - min(data(1, :)) + window;
x = data(2, :) - min(data(2, :)) + window;
initPos = round([y;x])';

imshow(cIM, []); hold on;
% range = getrangefromclass(cIM);
% waitbar(0);

for i = 1:numel(x)
    out(i).BWx = [];
    out(i).BWy = [];
end

% print('Par Loop start')
% warning('off','all')
% parfor i = 1:numel(x)
for i = 1:20
    imCut = cIM(initPos(i,1)-25:initPos(i,1)+25, initPos(i,2)-25:initPos(i,2)+25);
    scatter(initPos(i,2),initPos(i,1), 40, 'g', 'fill'); hold on;
% end
%     thresh = graythresh(imCut);
%     thresh_bw = thresh*range(2);
%     thresh_bw = thresh_bw - thresh_bw*0.5;
    
    thresh_bw = mean(imCut(:)) - std(double(imCut(:)))/10 ;
    if thresh_bw < 0
        thresh_bw = 0;
    end
    
    indS = find(imCut>quantile(imCut(:), 0.99));
    if numel(indS) > 0
        pick = round(numel(indS)/2);
        [ind1, ind2] = ind2sub([51,51], indS(pick));
        
        RGpos = [ind1 + initPos(i,1)-26, ind2 + initPos(i,2)-26];
        
        cIM(RGpos(1), RGpos(2));
        [~, BW] = regionGrowing(cIM, RGpos, thresh_bw, 15);
        [BWx,BWy] = ind2sub(size(cIM),find(BW>0));
        scatter(BWy, BWx, 20, 'b', 'fill')
        out(i).BWx = BWx +  min(data(1, :)) - window;
        out(i).BWy = BWy +  min(data(2, :)) - window;
    end
    
    % waitbar(i/numel(x))
    
    hold on;
    fprintf('%.0f/%.0f\n', i, numel(x));
%     print(strcat(num2str(i), '/', num2str(numel(x))));     
end
% warning('on','all')
save(save_p, 'out')  
dummy = 0;
end


