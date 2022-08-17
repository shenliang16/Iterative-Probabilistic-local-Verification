function [Crrsp, time, itr, InlierRatio, InlierRatio_inlier, InlierRatio_outlier] = IPLV_Rigid(X, Y, Sift_Ratio, omiga, Alpha, sigma, NeiNum, thr, itrMax, flag, Crr_GT) 

%% flag:   -1  : neighborhood_Test
%% flag:   0  : KNN  
%% flag:   1  : CDR   1 iteration
%% flag:   2  : CDR   2 iteration

N=size(X,1);  j=0;   D=2;  M=1;             InlierRatio=[]; InlierRatio_inlier=[]; InlierRatio_outlier=[];

%% parameters setting  
if ~exist('NeiNum','var');      NeiNum = [8,8];   end
if ~exist('omiga','var');       omiga=0.2;           end
if ~exist('Alpha','var');       Alpha=0.3;          end
if ~exist('sigma','var');       sigma=0.02;          end
if ~exist('dis_ratios','var');  dis_ratios=ones(N,1);          end
if ~exist('thr','var');         thr=0.5;          end

if length(NeiNum)==1;      NeiNum = [NeiNum,NeiNum];   end
K1 = NeiNum(1);     K2 = NeiNum(2); 
omiga = omiga*ones(N,1);


%%  iteration 1
tic; 
X1 = X(:,1:2);   Y1 = Y(:,1:2);  
[X1,Y1,normal] = AffineCPD_normalize(X1,Y1);


%% Scaling and Rotation
Scaling = Y(:,3)./X(:,3);               Scaling = Scaling*normal.xscale/normal.yscale; 
Theta = Y(:,4) - X(:,4);
% Rotation = [cos(Theta), sin(Theta); -sin(Theta), cos(Theta)];


%% Rigid Matching
% K1 = round(K1*Alpha);
% Alpha = [];
% K1 = 1;
Pcoe = 1;
ff = exp(-Pcoe*Sift_Ratio);


if flag<=0
%% KNN 
    NumNei=K1+1; 
    Data=X1(:, 1:2)'; kdtreeX = vl_kdtreebuild(Data);
    [neighborX, ~] = vl_kdtreequery(kdtreeX, Data, Data, 'NumNeighbors', NumNei) ;
    if flag == -1
        % 邻域质量  KNN_X
        [InlierRatio(1), InlierRatio_inlier(1), InlierRatio_outlier(1), ~] = NeigInlierRatio3(neighborX,  Crr_GT);
        NumNei=K1+1; 
        Data=Y1(:, 1:2)'; kdtreeX = vl_kdtreebuild(Data);
        [neighborX, ~] = vl_kdtreequery(kdtreeX, Data, Data, 'NumNeighbors', NumNei) ;
        % 邻域质量  KNN_Y
        [InlierRatio(2), InlierRatio_inlier(2), InlierRatio_outlier(2), ~] = NeigInlierRatio3(neighborX,  Crr_GT);
    end
end

if flag~=0
%% CDR
    Kmax = 4*K1;%Kmax=round(max(min(100,(N-1)/2),K1));%ff=(ff-min(ff))./(max(ff)-min(ff));
    Xt=[X(:, 1:2)']; kdtreeX = vl_kdtreebuild(Xt); 
    [neighborX_temp1, DisX1] = vl_kdtreequery(kdtreeX, Xt, Xt, 'NumNeighbors', Kmax+1);
    Xt=[Y(:, 1:2)']; kdtreeX = vl_kdtreebuild(Xt); 
    [neighborX_temp2, DisX2] = vl_kdtreequery(kdtreeX, Xt, Xt, 'NumNeighbors', Kmax+1);
    neighborX_temp = [neighborX_temp1; neighborX_temp2(2:end, :)];
    DisX = [DisX1; DisX2(2:end, :)];
    f = 1-Sift_Ratio; ff=min(f,0.2); 
    Dis = DisX(2:end, :);
    Sp_nei_Kmax=neighborX_temp(2:end,:); % 保留最近邻域（局部区域）信息，一边更新中使用
    Sp_nei=Sp_nei_Kmax;  %Dis=max(Dis,Dis(4*K1));%Dis=max(Dis,ReluPara);  %Dis = sigmf(Dis,[1/10 ReluPara]);%
    f_local=ff(Sp_nei)./(Dis);   %局部置信度排名
    [~,ind]=sort(f_local,'descend');
    for i=1:N; Sp_nei(1:K1,i)=Sp_nei(ind(1:K1,i),i); end
    Sp_nei = Sp_nei(1:K1,:);
    neighborX = [neighborX_temp(1, :); Sp_nei];
    % 邻域质量
    if flag == -1
        [InlierRatio(3), InlierRatio_inlier(3), InlierRatio_outlier(3), ~] = NeigInlierRatio3(neighborX,  Crr_GT);
    end
end


%% First iteration  仅CDR
if flag == -1 || flag==2
    %% 预计算
    DisX2 = (X1(neighborX(2:end,:), 1:2) - X1(repmat(neighborX(1,:),[K1,1]), 1:2))';     DisX = reshape(DisX2,2,K1,[]);% 中心化
    DisY2 = (Y1(neighborX(2:end,:), 1:2) - Y1(repmat(neighborX(1,:),[K1,1]), 1:2))';     DisY = reshape(DisY2,2,K1,[]);% 中心化     
    Pu = 10;
    c=zeros(N,1);   c2=zeros(K1,N);  sontemp1=c2;    ErrorSum=zeros(N,K1);    distL =zeros(K1,D,N);
    for i = 1 : N
        % local 
        Rotation = [cos(Theta(i)), sin(Theta(i)); -sin(Theta(i)), cos(Theta(i))];   % Rotation
        Error = (DisY(:,:,i) - Scaling(i) .* (Rotation * DisX(:,:,i)))';                      % Error
        [ErrorSort, ind] = sort(sum(Error.^2,2));
        ErrorSum(i,:) = ErrorSort(1:K1);
        distL(:,:,i) = Error(ind(1:K1),:).^2;    % 计算了两次平方，可优化
    end
    %%
    if size(omiga)==1
       omiga = omiga*ones(N,1);
    end
    for i = 1 : N
        [c(i,1),c2(:,i),sontemp1(:,i)] = BayesEstimatorL(distL(:,:,i),Pu,omiga(i),Alpha,sigma);
    end
    Kmax = 6*K1;%Kmax=round(max(min(100,(N-1)/2),K1));%ff=(ff-min(ff))./(max(ff)-min(ff));
    Xt=[X(:, 1:2)']; kdtreeX = vl_kdtreebuild(Xt); 
    [neighborX_temp1, DisX1] = vl_kdtreequery(kdtreeX, Xt, Xt, 'NumNeighbors', Kmax+1);
    Xt=[Y(:, 1:2)']; kdtreeX = vl_kdtreebuild(Xt); 
    [neighborX_temp2, DisX2] = vl_kdtreequery(kdtreeX, Xt, Xt, 'NumNeighbors', Kmax+1);
    neighborX_temp = [neighborX_temp1; neighborX_temp2(2:end, :)];
    DisX = [DisX1; DisX2(2:end, :)];
    f = c; ff=min(f,0.5); 
    Dis = DisX(2:end, :);
    Sp_nei_Kmax=neighborX_temp(2:end,:); % 保留最近邻域（局部区域）信息，一边更新中使用
    Sp_nei=Sp_nei_Kmax;  %Dis=max(Dis,Dis(4*K1));%Dis=max(Dis,ReluPara);  %Dis = sigmf(Dis,[1/10 ReluPara]);%
    f_local=ff(Sp_nei)./(Dis);   %局部置信度排名
    [~,ind]=sort(f_local,'descend');
    for i=1:N; Sp_nei(1:K1,i)=Sp_nei(ind(1:K1,i),i); end
    Sp_nei = Sp_nei(1:K1,:);
    neighborX = [neighborX_temp(1, :); Sp_nei];
    % 邻域质量
    if flag == -1
        [InlierRatio(4), InlierRatio_inlier(4), InlierRatio_outlier(4), ~] = NeigInlierRatio3(neighborX,  Crr_GT);  
    end
end



%% 预计算
DisX2 = (X1(neighborX(2:end,:), 1:2) - X1(repmat(neighborX(1,:),[K1,1]), 1:2))';     DisX = reshape(DisX2,2,K1,[]);% 中心化
DisY2 = (Y1(neighborX(2:end,:), 1:2) - Y1(repmat(neighborX(1,:),[K1,1]), 1:2))';     DisY = reshape(DisY2,2,K1,[]);% 中心化       
Pu = 10;
c=zeros(N,1);   c2=zeros(K1,N);  sontemp1=c2;    ErrorSum=zeros(N,K1);    distL =zeros(K1,D,N);
for i = 1 : N
    % local 
    Rotation = [cos(Theta(i)), sin(Theta(i)); -sin(Theta(i)), cos(Theta(i))];   % Rotation
    Error = (DisY(:,:,i) - Scaling(i) .* (Rotation * DisX(:,:,i)))';                      % Error
    [ErrorSort, ind] = sort(sum(Error.^2,2));
    ErrorSum(i,:) = ErrorSort(1:K1);
    distL(:,:,i) = Error(ind(1:K1),:).^2;    % 计算了两次平方，可优化
end

%%
for itr = 1:itrMax
    
    if itr>1
        Np = sum(c);        
        DPsum = c'*sum(c2'.*ErrorSum,2); 
        Psum = c'*sum(c2)';

        Q(itr) = DPsum/(2*sigma) + Np*log(sigma) - sum(c.*log(omiga)) - sum((1-c).*log(1-omiga));

        omiga = Np/N;                                         %omiga=min(omiga,0.8);
%         omiga = omiga0 * (Np/Sum_omiga);                       %omiga=min(omiga,0.8);
        Alpha=Psum/(K1*Np);     %Alpha=sum(c2)/K1;
        sigma = DPsum/(2*Psum);               
        Alpha1(itr)=mean(Alpha);    sigma1(itr)=sigma;   omiga1(itr)=omiga;    %     fprintf(sprintf('omiga:%.3f;    Alpha:%.3f;    sigma:%.3f;    S:%.3f;\n', omiga,  Alpha,  sigma,  S));

        if abs((Q(itr)-Q(itr-1))/Q(itr))<1e-1  || sigma<1e-6
            break;
        end
    end
    
    if size(omiga)==1
       omiga = omiga*ones(N,1);
    end
    if size(Alpha)==1
       Alpha = Alpha*ones(N,1);
    end
    
    for i = 1 : N
        [c(i,1),c2(:,i),sontemp1(:,i)] = BayesEstimatorL(distL(:,:,i),Pu,omiga(i),Alpha(i),sigma);
    end
end

% figure; plot(sigma1);   figure; plot(Alpha1);   figure; plot(Q);   figure; plot(sort(c))


ind_crr= find(c>thr(1));   %itr=j;

%% 输出
time=toc;
% [InlierRatio,~]=NeigInlierRatio(neighbortemp');
Crrsp = ind_crr';
if isempty(Crrsp)
    Crrsp=[1:8];
end
% Crrsp = [Crrsp;Crrsp]';