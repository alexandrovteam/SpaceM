function [X,vLabels]=Generate_MNIST(NofPts, mOpts) 

%
% function [X,vLabels]=Generate_MNIST(NofPts, mOpts) 
% 
% Generate_MNIST gets handwritten digits from the MNIST data set.
%
% IN: 
%    NofPts: a number or a vector, the number of points to extract. The vector case, the number of points from each digits. See QueryPts in the Opts.
%    [mOpts]: structure containing the following fields:
%             [Sampling]: 'FirstN' returns the first NofPts, 'RandN' returns the NofPts points after a random permutation. Default = RandN
%             [QueryDigits]: if it is empty, it doesn't distinguish digits. If it is a vector cosisting of digits, only those digits are sampled. 
%                         ex) [0, 3, 5]. In this case, the NofPts should be matched. Otherwise, it samples NofPts(1) points from each digit.
%                         Default =[].
%             [ReturnForm]: if 'vector', it returns a point as a 28^2 x 1 vector. if 'matrix', it returns as a 28x28 matrix. Defalut = 'vector' 
%
% OUT: a double sum(NofPts)x28^2 matrix if 'vector' is the ReturnForm, a double sum(NofPts)x28x28 3-dim array if 'matrix'. 
%
% Example:  X=Generate_MNIST([500, 500, 500], struct('Sampling', 'RandN', 'QueryDigits', [1, 6, 9], 'ReturnForm', 'matrix'));
%           X=Generate_MNIST(1000);
%
% SC:
%    YM: 8/20/2008
%    MM: 12/2/2008
%

%
% (c) Copyright Duke University, 2008
% Mauro Maggioni
%

% Setup parameters
if nargin < 2
    mOpts = [];
end

if ~isfield(mOpts, 'Sampling')
    mOpts.Sampling='RandN';
end
if ~isfield(mOpts, 'QueryDigits')
    mOpts.QueryDigits=[];
end
if ~isfield(mOpts, 'ReturnForm')
    mOpts.ReturnForm='vector';
end


% Loading data
load TrainImages;
load TrainImageLabels;


% Getting digits
if isempty(mOpts.QueryDigits)
    if NofPts > length(Labels)
        fprintf('\n The number of points is larger than the number in the data set. It will return the full data set.');
        NofPts=length(Labels);
    end        
    if  strcmpi(mOpts.Sampling, 'FirstN')
        X=TrainImages(1:NofPts, :, :);       
        vLabels = Labels(1:NofPts)';
    else 
        idxs=randperm(length(Labels));
        idxs=idxs(1:NofPts);
        X=TrainImages(idxs, :, :);
        vLabels = Labels(idxs)';
    end
else 
    QueryDigits=mOpts.QueryDigits;
    NofQPts=length(QueryDigits);
    if NofQPts > 1
        if length(NofPts)==1
            %fprintf(['\n The ', num2str(NofPts), ' will be sampled from each digit.']);
            NofPts=NofPts*ones(1, NofQPts);           
        elseif length(NofPts)~=NofQPts
            %fprintf(['\n The size of points and query digits are different. ', num2str(NofPts(1)), ' points will be sampled from each digit.']);
            NofPts=NofPts(1)*ones(1, NofQPts);
        end
    end
    
    % seperating digits
    clabels=cell(1,10);
    for i=1:10
        clabels{i}=find(Labels==i-1); % note 0 is label1, 1 is label2,.....
    end
    LengthL=zeros(10,1);
    for i=1:10
        LengthL(i)=length(clabels{i});
    end   

    if NofQPts==1  
        if NofPts > LengthL(QueryDigits+1)
            fprintf(['\n The number of points is larger than the number in the digit', num2str(QueryDigits), '. It will return the full data set.']);
            NofPts=length(Labels);
        end
        if  strcmpi(mOpts.Sampling, 'FirstN')
            X=TrainImages(clabels{QueryDigits+1}(1:NofPts), :, :);  
            vLabels=Labels(clabels{QueryDigits+1}(1:NofPts))';
        else 
            idxs=randperm(length(clabels{QueryDigits+1}));
            idxs=idxs(1:min([NofPts,length(clabels{QueryDigits+1})]));
            X=TrainImages(clabels{QueryDigits+1}(idxs), :, :);
            vLabels=Labels(clabels{QueryDigits+1}(idxs))';
        end        
    else  
        for i=1:NofQPts
            if NofPts(i) > LengthL(QueryDigits(i)+1)
                fprintf(['\n The number of points is larger than the number in the digit', num2str(QueryDigits(i)), '. It will return the full data set.']);
                NofPts(i)=LengthL(QueryDigits(i)+1);
            end
        end
        X=zeros(sum(NofPts), size(TrainImages, 2), size(TrainImages,3));
        vLabels = zeros(sum(NofPts),1);
        sumNofPts=0;
        for i=1:NofQPts         
            if  strcmpi(mOpts.Sampling, 'FirstN')
                X(sumNofPts+1:sumNofPts+NofPts(i), :, :) = TrainImages(clabels{QueryDigits(i)+1}(1:NofPts(i)), :, :);  
                vLabels(sumNofPts+1:sumNofPts+NofPts(i)) = Labels(clabels{QueryDigits(i)+1}(1:NofPts(i)));
            else 
                idxs=randperm(length(clabels{QueryDigits(i)+1}));
                idxs=idxs(1:NofPts(i));
                X(sumNofPts+1:sumNofPts+NofPts(i), :, :) = TrainImages(clabels{QueryDigits(i)+1}(idxs), :, :);
                vLabels(sumNofPts+1:sumNofPts+NofPts(i)) = Labels(clabels{QueryDigits(i)+1}(idxs));
            end
            sumNofPts=sumNofPts+NofPts(i);
        end
    end 
end

% if  length(mOpts.QueryDigits)<=1
%     NofQPts=1;
% end
% sumNofPts=0;
% for i=1:NofQPts
%     figure;
%     imshow(TileBlocks(X, sumNofPts+1:sumNofPts+NofPts(i), min(NofPts(i), 20)));
%     sumNofPts=sumNofPts+NofPts(i);
% end

if strcmpi(mOpts.ReturnForm, 'vector')
    X=reshape(X, [size(X, 1), size(X,2)*size(X,3), 1]);
end
X=double(X);


return;














