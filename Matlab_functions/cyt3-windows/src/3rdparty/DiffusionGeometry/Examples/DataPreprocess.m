function Data=DataPreprocess(Data, POpts)

%
% function Data=DataPreprocess(Data, POpts)
%
% preprocess data set using convolution, regularization, etc
%
% input: 
%   Data : double L by M array, L numbers of M dim vector
%          or  double L by M by N array, L numbers of M by N matrix 
%   POpts: process options
%       [Conv]: if true, convolution using a kernel with entry 1
%           [COnvKSize]: size of convolution kernel, default=3
%           [ConvKSum]: if 1, kernel is normalized as sum 1
%           [ConvN]: the number of convolution, default=1
%       [LapLeg]: if true, regularize data using Laplacian with a data
%                 fitting term. 
%                 Rmk: Do not use Conv and LapLeg both. Roles are similar
%           [LapLamda]: parameter for data fitting term, defalut = 1
%           [LapIterN]: iteration number, default = 3                     
%       [MeanZero]: to get mean zero data
%       [UnitBall]: to locate all data in the unit ball
%       [UnitSphere]: to locate all nonzero data on the unit sphere
%           [UnitSphereToler]: if the norm is smaller than this number, 
%                              it is located at the origin 0. default = 0
%
% output:
%   Data: preprocessed data, L by M array or L by M by N array
% 
% use:
%   LapReg1D, LapReg2D
%
% SC: YM     6/19/2008 
%


if nargin==1 || isempty(POpts)
    fprintf('No processing options are given. So no preprocess.\n');
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Collecting options....

if isfield(POpts, 'Conv') && POpts.Conv
    Conv=true;
    if isfield(POpts, 'ConvKSize')
        ConvKSize=POPts.ConvKSize;
    else 
        ConvKSize=3;
        fprintf('Convolution kernel size is not given, 3 is assgined.\n');
    end 
    if isfield(POpts, 'ConvKSum') 
        ConvKSum =Opts.ConvKSum;
    else 
        ConvKSum =ConvKSize^2;
    end
    if isfield(POpts, 'ConvN')
        ConvN=POtps.ConvN;
    else
        ConvN=1;
    end
else
    Conv=false;
end

if isfield(POpts, 'LapReg') && POpts.LapReg
    LapReg=true;
    if isfield(POpts, 'LapLambda')
        Lambda=POPts.LapLambda;
    else 
        Lambda=1;
        fprintf('Lambda for Laplacian regularization is not given, 1 is assigned.\n');
    end 
    if isfield(POpts, 'LapIterN')
        LapIterN=POpts.IterN;
    else
        LapIterN=3;
        fprintf('Iteration number for Laplacian regularization is not given, 3 is assigned.\n');
    end
else
    LapReg=false;
end

if isfield(POpts, 'MeanZero') && POpts.MeanZero
    MeanZero=true;
else
    MeanZero=false;
end

if isfield(POpts, 'UnitBall') && POpts.UnitBall
    UnitBall=true;
else
    UnitBall=false;
end

if isfield(POpts, 'UnitSphere') && POpts.UnitSphere
    UnitSphere=true;
    if isfield(POpts, 'UnitSphereToler') 
        UnitSphereToler=POpts.UnitSphereToler;
    else
        UnitSphereToler = 0;
    end
else
    UnitSphere=false;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preprocessing

if length(size(Data))==2  % vector case
    [L, M]=size(Data);
    
    if Conv
        ConvKer=ConvKSum/ConvKSize*ones(1, ConvKSize);
        for l=1:L,
            for li = 1:ConvN
                Data(l,:) = convn(Data(l,:), ConvKer, 'same');
            end
        end
    end

    
    if LapReg
        Data = LapReg1D(Data, Lambda, LapIterN);
    end
    
    
    if MeanZero    
        Mean=mean(Data, 2);
        Data = Data-repmat(Mean,[1, M]);
    end
    
    
    if UnitBall
        MaxL2norm=sqrt(max(sum(Data.^2, 2)));
        if MaxL2norm>0
            Data=Data/MaxL2norm;
        else
            error('All data are zero.\n');
        end
    end
    
    
    if UnitSphere
       L2norm=sqrt(sum(Data.^2, 2));
       Nzidxs=find(L2norm>UnitSphereToler);
       Data(Nzidxs, :)=Data(Nzidxs,:)./repmat(L2norm(Nzidxs), [1, M]);
    end
    
    
else  % matrix case  
    [L, M, N]=size(Data);
    
    % convolution
    if Conv
        ConvKer=ConvKSum/ConvKSize^2*ones(ConvKSize);
        for l=1:L,
            X = squeeze(Data(l,:,:));
            for li = 1:ConvN
                Data(l,:,:) = conv2(X, ConvKer, 'same');
            end;
        end;
    end

    
    %Laplacian Regularization
    if LapReg
        for l = 1:L,
            Data(l,:,:) = LapReg2D(squeeze(Data(l,:,:)), Lambda, LapIterN);
        end
    end

    
    if MeanZero       
        for l=1:L
            Mean = sum(sum(Data(l,:,:),3),2)/(M*N);
            Data(l,:,:) = Data(l,:,:)-Mean;
        end
    end

    
    if UnitBall
        MaxL2norm=sqrt(max(sum(sum(Data(:,:,:).^2,3),2)));
        if MaxL2norm>0
            Data=Data/MaxL2norm;
        else
            error('All data are zero.\n');
        end
    end
    
    
    if UnitSphere
       L2norm=sqrt(sum(sum(Data(:,:,:).^2,3),2));
       Nzidxs=find(L2norm>UnitSphereToler);
       NofNz=length(Nzidxs);
       for i=1:NofNz
           Data(Nzidxs(i),:,:)=Data(Nzidxs(i), :, :)/L2norm(Nzidxs(i));
       end
    end   
        
end

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function U=LapReg1D(U0, lambda, stop)

%
% function U=LapReg1D(U0, lambda, stop)
%
% Tychonoff regularization using grad^2
% for 1-D with Neumann condition
% Gauss-Jacobi iterative method
%
% input: 
%   U: a row vector or matrix, denoising each row, size M by N
%   stop: the number of iteration
%   lambda: weight on the fidelity term
%
% ouput: 
%    U_new: denoised row vector or matrix, size M by N
% Date: Mar 5, 2008
%

[M,N]=size(U0); 

%Preprocessing: removing Inf and NaN
badPts=find((isinf(U0))|(isnan(U0)));
U0(badPts)=0;
Lambda=lambda*ones(M,N);
Lambda(badPts)=0;


U=U0;

for i=1:stop,
    Uext=[U(:,1), U, U(:, N)]; % Neumann extension, M by N+2 matrix
    U=(Uext(:,1:N)+Uext(:,3:N+2)+Lambda.*U0)./(2+Lambda);
end

return;



function U=LapReg2D(U0, lambda, stop)

%
% function U=LapReg2D(U0, lambda, stop)
%
% Regurizing or smoothing  using
% int |grad U|^2 + lambda|U-U0|^2 
% If lambda = 0, it is the heat defusion 
% Boundary condition: Newmann cond
% grid is unform and rectangular
%
% Input:
%   U0: given image (double, [0, 1])
%   lambda: fidelity const
%   stop: iteration number
%
% Output:
%   U: regularized image (double, [0, 1])
%

[m, n]=size(U0);

% Intializing
U=U0;

for t=1:stop
    Uext=[U(:,1), U, U(:, n)]; %Neumann extension in x direction, M by N+2 matrix
    Uext=[Uext(1, :); Uext ;Uext(m,:)]; %Neumann extension in y direction, M+2 by N+2 matrix
    % By extension, index (i, j) in Uext is (i+1, j+1)
    U=(Uext(2:m+1,1:n)+Uext(2:m+1, 3:n+2)+Uext(1:m, 2:n+1)+Uext(3:m+2, 2:n+1)+lambda*U0)./(4+lambda);    
end

return;


















