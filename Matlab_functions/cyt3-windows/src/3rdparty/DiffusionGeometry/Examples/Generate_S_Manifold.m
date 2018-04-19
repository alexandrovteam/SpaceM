function cX=Generate_S_Manifold(NofPts, Dim, Opts)

%
% function cX=Generate_S_Manifold(NofPts, Dim, Opts)
%
% Generate_S_Manifold generates a S shape manifold.
% 
% IN:
%    NofPts     : the number of points in the manifold generated
%    [Dim]      : the dimension of the manifold, if Dim=1, a curve, if Dim=2, a surface. default = 2
%    [Opts]     : structure containing the following fields:
%                   [PtsType] : 'mesh': a meshgrid type 
%                               'rand': a random uniform sampling
%                               default = rand
%
% OUT:
%     cX: NofPtsxDim array, if PtsType is mesh, not exactly. 
%
% Example: X = Generate_S_Manifold(1000, 2, struct('PtsType', 'mesh'));
%          X = Generate_S_Manifold(1000);
%
% SC:
%    YM: 8/19/2008
%


% Setup parameters
if nargin < 2
    Dim = 2;
end

if nargin < 3
   Opts=[];
end

if Dim > 2
    fprintf('Dim > 2, Dim is modified to 2.\n'); 
    Dim = 2;
end

if ~isfield(Opts, 'PtsType')
    Opts.PtsType = 'rand';
end

if strcmpi(Opts.PtsType,'mesh') && Dim ==2
    fprintf('The mesh type case, it may not return the exact number of points.\n');
end

% Generate data
if Dim == 1    % curve case
    if strcmpi(Opts.PtsType, 'mesh') % mesh type sampling
        if mod(NofPts, 2)
            HalfNofPts=floor(NofPts/2)+1;
            dx1=linspace(0, 1, HalfNofPts);
            dx2=linspace(0, 1, HalfNofPts)+1/(HalfNofPts-1);
            dx2(HalfNofPts)=[];
        else           
            HalfNofPts=NofPts/2;
            dx1=linspace(0, 1, HalfNofPts);
            dx2=linspace(0, 1, HalfNofPts)+1/(HalfNofPts-1);
        end
    else    % random sampling from uniform dist
        dx=rand(1, NofPts);
        if mod(NofPts, 2)
           HalfNofPts=floor(NofPts/2)+1;
        else
           HalfNofPts=NofPts/2; 
        end
        dx1=dx(1:HalfNofPts);
        dx2=dx(HalfNofPts+1:NofPts);
    end
    dx1=pi/2+3*pi/2*dx1;
    dx2=pi-3*pi/2*dx2;
    X1=[cos(dx1')-1, sin(dx1')]/2;
    X2=[cos(dx2')+1, sin(dx2')]/2;
    cX=[X1;X2];
%     size(cX, 1)
%     figure; plot(cX(:, 1), cX(:, 2), '.');axis equal

elseif Dim == 2 % surface case
    if strcmpi(Opts.PtsType, 'mesh') % meshtype sampling
%         GridN=round(sqrt(NofPts/2));        
%         dx1=linspace(0,1, GridN); 
%         dx2=linspace(0,1, GridN)+1/(GridN-1);
%         dz=linspace(0,1, GridN);
%         dx1=pi/2+3*pi/2*dx1;
%         dx2=pi-3*pi/2*dx2;
%         dxdy1=[cos(dx1)-1; sin(dx1)]/2;
%         dxdy2=[cos(dx2)+1; sin(dx2)]/2;
%         dxdy=[dxdy1,dxdy2];
%         [meshdx, meshdz]=meshgrid(dxdy(1, :), dz);
%         [meshdy, meshdz]=meshgrid(dxdy(2, :), dz);
%         meshdx=reshape(meshdx, 1, 2*GridN^2);
%         meshdy=reshape(meshdy, 1, 2*GridN^2);
%         meshdz=reshape(meshdz, 1, 2*GridN^2);
%         cX=[meshdx;meshdy;meshdz]';    

        GridN=round(sqrt(NofPts/9)); 
           if mod(9*GridN, 2)
            Mid=floor(9*GridN/2)+1;
            dx1=linspace(0, 1, Mid);
            dx2=linspace(0, 1, Mid)+1/(Mid-1);
            dx2(Mid)=[];
        else           
            Mid=9*GridN/2;
            dx1=linspace(0, 1, Mid);
            dx2=linspace(0, 1, Mid)+1/(Mid-1);
        end
        dz=linspace(0,1, GridN);
        dx1=pi/2+3*pi/2*dx1;
        dx2=pi-3*pi/2*dx2;
        dxdy1=[cos(dx1)-1; sin(dx1)]/2;
        dxdy2=[cos(dx2)+1; sin(dx2)]/2;
        dxdy=[dxdy1,dxdy2];
        [meshdx, meshdz]=meshgrid(dxdy(1, :), dz);
        [meshdy, meshdz]=meshgrid(dxdy(2, :), dz);
        meshdx=reshape(meshdx, 1, 9*GridN^2);
        meshdy=reshape(meshdy, 1, 9*GridN^2);
        meshdz=reshape(meshdz, 1, 9*GridN^2);
        cX=[meshdx;meshdy;meshdz]';        


    else   % random sampling from uniform dist
        dx=rand(1, NofPts);
        X3=rand(NofPts, 1);       
        if mod(NofPts, 2)
           HalfNofPts=floor(NofPts/2)+1;
        else
           HalfNofPts=NofPts/2; 
        end
        dx1=dx(1:HalfNofPts);
        dx2=dx(HalfNofPts+1:NofPts);      
        dx1=pi/2+3*pi/2*dx1;
        dx2=pi-3*pi/2*dx2;
        X1=[cos(dx1')-1, sin(dx1')]/2;
        X2=[cos(dx2')+1, sin(dx2')]/2;
        cX=[[X1; X2], X3];   
    end 
%     size(cX, 1)
%     figure; plot3(cX(:, 1), cX(:, 2), cX(:, 3), '.');axis equal
end
    
return;

