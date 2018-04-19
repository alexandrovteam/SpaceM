function [G,Coords]=FilterGraph(A,par,level,NN)

%
% function G=filtergraph(A,par,level,NN);
%
% IN:
%   A         : an image of dimensions H by W
%   par       : one of:
%       'ptch'  a patch graph, level determines size of patch
%       'rand'  white noise graph, NN determines size of white noise filter,
%                level determines number of instances
%       'haar'  haar functions.  level determines number of scales
%       'vari'  variances in a patch.  level determines size of path, NN
%               the number of moments; so level=5 and NN=3 gives a graph with
%               the first coordinate as the 5x5 averages, the second coordinate
%               the variances over the patches, and the third and fourth coordinate
%               the third and fourth moments, etc.
%       'unif'  same as rand but with uniform noise patches.  not normalized.
%       'loco'  dct filters.  level determines the size of the patch;  same
%               as patch up to isometry.
% OUT:  
%   G, a (HxW,d) matrix, where d depends on the paramters.
%
% note: uses image processing toolbox

% SC:
%   MM:     3/4/2011
%

G       = [];
Coords  = [];

% N is the (odd) side length of the cube

%if size(A,2)==3
%    A=cubify(A(:,3),A,H,W);
%end

[H,W] = size(A);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A=double(A);
if par=='haar'
    G=zeros(H*W,level*3+1);
    for k=level:-1:1,
        d=ones(2^(k-1));

        Hm=[d d; -d -d];

        Ha{k}(1).H=(Hm)/norm(Hm,'fro');

        Hm=[d -d; d -d];

        Ha{k}(2).H=(Hm)/norm(Hm,'fro');

        Hm=[d -d; -d d];

        Ha{k}(3).H=(Hm)/norm(Hm,'fro');

    end
    Hm=ones(2^level);
    Ha{level+1}(1).H=(Hm)/norm(Hm,'fro');


    for l=1:level
        for s=1:3

            %Q=imfilter(A,Ha{l}(s).H);
            Q=conv2(A,Ha{l}(s).H,'same');

            %B=zeros(H*W,1);
            % for k=1:H
            % for j=1:W
            %     a=Q(k,j);
            %     B((j-1)*W+k,1)=a;
            % end
            % end

            % G(:,(l-1)*3+s)=B;
            G(:,(l-1)*3+s)=reshape(Q,H*W,1);

        end
    end
    Q=conv2(A,Ha{level+1}(1).H,'same');
    %Q=imfilter(A,Ha{level+1}(1).H);

    %B=zeros(H*W,1);
    %for k=1:H
    %for j=1:W
    %    a=Q(k,j);
    %    B((j-1)*W+k)=a;
    %end
    %end
    %G(:,s*level+1)=B;
    G(:,s*level+1)=reshape(Q,H*W,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

elseif par=='edge'
    %level=2;
    %NN=81;
    for k=1:level

        Ha{k}=zeros(NN);
        for s=1:NN
            for t=1:NN
                if s~=ceil(NN/2)||t~=ceil(NN/2)
                    v1=[(s-ceil(NN/2)) (t-ceil(NN/2))];
                    v1=v1/norm(v1);
                    %v2=[cos(2*pi*(k-1)/level); sin(2*pi*(k-1)/level)];
                    v2=[cos(pi*(k-1)/level); sin(pi*(k-1)/level)];
                    a=abs(v1*v2)^2;
                    b=v1(1)*v2(2)-v1(2)*v2(1);
                    if abs(b)>.001
                        b=b/abs(b);
                    else
                        b=0;
                    end

                else
                    b=0;
                    a=1;
                end
                rrr(s,t,k)=b;
                Ha{k}(s,t)=erf(4*b*(1-a));
            end
        end
    end

    for l=1:level

        Q=imfilter(A,Ha{l},'circular');

        %B=zeros(H*W,1);
        % for k=1:H
        % for j=1:W
        %     a=Q(k,j);
        %     B((j-1)*W+k,1)=a;
        % end
        % end


        %G(:,l)=B;
        G(:,l)=reshape(Q,H*W,1);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif par=='loco'
    order=[1 3  6  11 18 27 38 51 66;
        2  4  8  13 20 29 40 53 68;
        5  7  9  15 22 31 42 55 70;
        10 12 14 16 24 33 44 57 72;
        17 19 21 23 25 35 46 59 74;
        26 28 30 32 34 36 48 61 76;
        37 39 41 43 45 47 49 63 78;
        50 52 54 56 58 60 62 64 80;
        65 67 69 71 73 75 77 79 81;];
    if level>9

        for j=1:level
            for k=1:level
                s=zeros(level);
                s(j,k)=1;
                Ha{(j-1)*level+k}=idct2(s)/norm(idct2(s),'fro');
            end
        end

    else
        order=order(1:level,1:level);

        for j=1:level^2
            [ii jj]=find(order==j);
            s=zeros(level);
            s(ii,jj)=1;
            Ha{j}=idct2(s)/norm(idct2(s),'fro');

        end
    end



    for l=1:level^2

        Q=imfilter(A,Ha{l});

        %B=zeros(H*W,1);
        % for k=1:H
        % for j=1:W
        %     a=Q(k,j);
        %     B((j-1)*W+k,1)=a;
        % end
        % end


        %G(:,l)=B;
        G(:,l)=reshape(Q,H*W,1);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

elseif par=='ptch'
    G=zeros(H*W,level^2);
    % Prepare the filters, which are just the elementary matrices
    for j=1:level
        for k=1:level
            s=zeros(level);
            s(j,k)=1;
            Ha{(k-1)*level+j}=s;
            %Ha{(k-1)*level+j}=imfilter(Ha{k},ones(3),'circular');
        end
    end

    for l=1:level^2
        Q=imfilter(A,Ha{l});
        G(:,l)=reshape(Q,H*W,1);        
        %B=zeros(H*W,1);
        % for k=1:H
        % for j=1:W
        %     a=Q(k,j);
        %     B((j-1)*W+k,1)=a;
        % end
        % end

        %G(:,l)=B;
    end;   
    [tmp1,tmp2]=meshgrid(1:H,1:W);
    Coords(:,1)=tmp1(:);
    Coords(:,2)=tmp2(:);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif par=='rand'
    if nargin==5
        NN=5;
    end
    for k=1:level

        Ha{k}=randn(NN);
        Ha{k}=Ha{k}/norm(Ha{k},'fro');
        %Ha{k}=imfilter(Ha{k},ones(3),'circular');
    end

    for l=1:level

        Q=imfilter(A,Ha{l},'circular');

        %B=zeros(H*W,1);
        %for k=1:H
        %for j=1:W
        %    a=Q(k,j);
        %    B((j-1)*W+k,1)=a;
        %end
        %end

        G(:,l)=reshape(Q,H*W,1);
        %G(:,l)=B;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

elseif par=='unif'

    for k=1:level

        Ha{k}=rand(5);
        %Ha{k}=imfilter(Ha{k},ones(3),'circular');
    end

    for l=1:level

        Q=imfilter(A,Ha{l},'circular');

        %B=zeros(H*W,1);
        % for k=1:H
        % for j=1:W
        %     a=Q(k,j);
        %     B((j-1)*W+k,1)=a;
        % end
        % end

        G(:,l)=reshape(Q,H*W,1);
        %G(:,l)=B;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif par=='mran'

    for k=1:level
        for j=1:NN
            M=randn(2^k+1);
            M=M./sqrt(sum(sum(M.^2)));
            Ha{(k-1)*NN+j}=M;
        end
    end
    for l=1:NN*level

        Q=imfilter(A,Ha{l},'circular');

        % B=zeros(H*W,1);
        %  for k=1:H
        %  for j=1:W
        %      a=Q(k,j);
        %      B((j-1)*W+k,1)=a;
        %  end
        %  end
        G(:,l)=reshape(Q,H*W,1);

        % G(:,l)=B;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

elseif par=='mdge'
    %level=5;
    %NN=10;
    for k=1:level
        a=linspace(-1,1,2^k+2);
        b=linspace(-1,1,2^k+1);
        c1=diff(exp(-NN*a.^2));
        c1=c1/max(c1);
        c2=exp(-NN*b.^2);
        c2=c2./max(c2);

        Ha{3*k-2}=zeros(2^k+1);
        Ha{3*k-1}=zeros(2^k+1);
        Ha{3*k}=zeros(2^k+1);
        for s=1:2^k+1
            for t=1:2^k+1
                Ha{3*k-2}(s,t)=c1(s)*c2(t);
                Ha{3*k-1}(s,t)=c2(s)*c1(t);
                Ha{3*k}(s,t)=c1(s)*c1(t);
            end
        end
        if k==level
            for s=1:2^level+1
                for t=1:2^level+1

                    Ha{3*k+1}(s,t)=c2(s)*c2(t);
                end
            end
        end
    end
    A=double(A);

    for l=1:3*level+1

        Q=imfilter(A,Ha{l},'circular');
        Q=Q/norm(Ha{l},'fro');
        %B=zeros(H*W,1);
        % for k=1:H
        % for j=1:W
        %     a=Q(k,j);
        %     B((j-1)*W+k,1)=a;
        % end
        % end
        G(:,l)=reshape(Q,H*W,1);

        %G(:,l)=B;
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif par=='vari'
    if nargin<6
        NN=1;
    end
    G=zeros(H*W,NN+1);
    GG=filtergraph(A,H,W,'ptch',level);
    for k=1:H*W
        G(k,1)=sum(GG(k,:))/level^2;
        for s=1:NN
            G(k,s+1)=sum((GG(k,:)-sum(GG(k,:))/level^2).^(s+1));
        end
    end

    %     nI=A;
    %     for k=1:H
    %         for j=1:W
    %             for s=1:level
    %                 for t=1:level
    %                     tt(s,t)=A(mod(k+s-ceil(level/2)-1,H)+1,mod(j+t-ceil(level/2)-1,W)+1);
    %                 end
    %             end
    %             ss=sum(sum((tt-sum(sum(tt))).^2));
    %             nI(k,j)=ss;
    %         end
    %     end

    %G=reshape(nI,H*W,1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % elseif par=='dcro'
    %
    %     order=[1 3  6  11 18 27 38 51 66;
    %      2  4  8  13 20 29 40 53 68;
    %      5  7  9  15 22 31 42 55 70;
    %      10 12 14 16 24 33 44 57 72;
    %      17 19 21 23 25 35 46 59 74;
    %      26 28 30 32 34 36 48 61 76;
    %      37 39 41 43 45 47 49 63 78;
    %      50 52 54 56 58 60 62 64 80;
    %      65 67 69 71 73 75 77 79 81;]
    %
    %
    %
    %
    %
    %
    %      if level>9
    %
    %         for j=1:level
    %             for k=1:level
    %                 s=zeros(level);
    %                 s(j,k)=1;
    %                 Ha{(j-1)*level+k}=idct2(s)/norm(idct2(s),'fro');
    %             end
    %         end
    %
    %     else
    %         order=order(1:level,1:level);
    %
    %         for j=1:level^2
    %             [ii jj]=find(order==j);
    %             s=zeros(level);
    %             s(ii,jj)=1;
    %             Ha{j}=idct2(s)/norm(idct2(s),'fro');
    %
    %         end
    %     end
    %
    %
    %
    %     for l=1:level^2
    %
    %             Q=imfilter(A,Ha{l});
    %
    %             B=zeros(H*W,1);
    %              for k=1:H
    %              for j=1:W
    %                  a=Q(k,j);
    %                  B((j-1)*W+k,1)=a;
    %              end
    %              end
    %
    %
    %              G(:,l)=B;
    %     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif par=='disk'
    cu=cubebuild(2*level+1,2);
    z=zeros((2*level+1)^2,1);
    for t=1:(2*level+1)^2
        z(t)=norm(cu(t,:)-[.5 .5]);
    end
    idx=find(z<=.5);
    cu=cu(idx,:);
    M=zeros(length(idx));
    for s=1:length(idx)
        for t=1:length(idx)
            M(s,t)=exp(-50*sum((cu(s,:)-cu(t,:)).^2));
        end
    end
    d=sum(M,2);
    d=1./sqrt(s);
    D=diag(d);
    M=D*(M*D);
    [V,E]=eig(M);
    V=fliplr(V);
    for k=1:length(idx);
        z=zeros((2*level+1)^2,1);
        z(idx)=V(:,k);
        Ha{k}=reshape(z,2*level+1,2*level+1);
    end



    for l=1:length(idx)

        Q=imfilter(A,Ha{l});

        %B=zeros(H*W,1);
        % for k=1:H
        % for j=1:W
        %     a=Q(k,j);
        %     B((j-1)*W+k,1)=a;
        % end
        % end


        %G(:,l)=B;
        G(:,l)=reshape(Q,H*W,1);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif par=='rdge'
    ee=1;
    %offset=0;
    %level=2;
    NN=7;
    for k=1:level
        theta=pi*(k-1)/(2*level);
        rot=[cos(theta) sin(theta);-sin(theta) cos(theta)];
        Ha{k}=zeros(NN);
        for s=1:NN
            for t=1:NN
                v=[(s-ceil(NN/2)); (t-ceil(NN/2))];
                v=rot*v;
                %rrr(s,t,k)=b;
                Ha{k}(s,t)=v(1)*exp(-(v(1)^2+v(2)^2)/ee)/ee;
            end
        end
    end

    for l=1:level

        Q=imfilter(A,Ha{l},'circular');

        %B=zeros(H*W,1);
        % for k=1:H
        % for j=1:W
        %     a=Q(k,j);
        %     B((j-1)*W+k,1)=a;
        % end
        % end


        %G(:,l)=B;
        G(:,l)=reshape(Q,H*W,1);
    end



end