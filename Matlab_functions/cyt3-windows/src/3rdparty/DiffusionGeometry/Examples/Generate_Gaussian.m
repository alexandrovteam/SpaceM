function cX=Generate_Gaussian(NofPts, Dim, StdDev)

%
% cX=Generate_Gaussian(NofPts, Dim, StdDev)
%
% Generate_Gaussian generates a Gaussian distribution in R^Dim space.
% 
% IN:
%    NofPts     : scalar, the number of points to be sampled
%    Dim        : scalar, the dimension of Gaussian distribution
%    [StdDev]   : scalar or 1xDim vector, standard deviation. 
%                 if scalar, std dev is const along all dimensions. Default: 1 
%
% OUT:
%     cX: NofPts x Dim array 
%
% Example: X = Generate_Gaussian(1000, 5)
%          X = Generate_Gaussian(1000, 3, [1,0.5, 0.25])
%
% SC:
%    YM: 8/27/2008
%


if nargin < 3
    StdDev = 1;
end

if length(StdDev)>1 && Dim~=length(StdDev)
    error('The length of StdDev does not match with given Dim!');
end

if length(StdDev)==1
    cX=StdDev*randn(NofPts,Dim);
else
    cX=randn(NofPts,Dim);
    for i=1:Dim
        cX(:, i)=StdDev(i)*cX(:, i);
    end
end

% figure;
% plot3(cX(:, 1), cX(:, 2), cX(:, 3), '.');axis equal;

return;
