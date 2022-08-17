
% Input
%   X1：             Feature frames of Image 1,  N*4 matrix, [location, scaling, rotation]
%   X2：             Feature frames of Image 2,  N*4 matrix, [location, scaling, rotation]
%   Sift_ratio:      Descriptor distance ratio values in the ratio test (NNDR method)

k=8*ones(1,2);  omiga_1=0.3;  omiga_2=0.5; sigma2=1e-3;  itrMax = 3;  flag = 2;       
[Clrv, time, iter] = IPLV_Rigid(X1, X2, Sift_ratio, omiga_1, omiga_2, sigma2, k, 0.5, itrMax, flag);
                