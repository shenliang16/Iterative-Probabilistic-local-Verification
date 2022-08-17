function [InlierRatio, InlierRatio_inlier, InlierRatio_outlier, ...
            InlierNum_all, InlierNum_inlier, InlierNum_outlier, Knn_ind_GT]=NeigInlierRatio3(Knn_ind,  Crr_GT_4, K1)
%% È¥³ıÖØ¸´
    [Nei_sort, index_sort] = sort(Knn_ind);
    temp1 = diff(Nei_sort);
    for i = 1:size(temp1, 1)
        for j = 1:size(temp1, 2)
            if temp1(i, j) == 0
                Knn_ind(index_sort(i, j)) = 0;
            end
        end
    end

%% È¥Áã
    Crr_GT_4 = [Crr_GT_4; 0];  
    if exist('N_all', 'var')
        N_all = sum(Knn_ind~=0)';    
    else
        N_all = K1;
    end
    Knn_ind(Knn_ind ==0) = length(Crr_GT_4);

%% 
    Knn_ind_GT = Crr_GT_4(Knn_ind)';
    InlierNum_all = sum(Knn_ind_GT(:, 2:end), 2);
    InlierRatio_all = InlierNum_all ./ N_all;

%% outliers
    ind = ~Knn_ind_GT(:, 1);
    InlierNum_outlier = InlierNum_all(ind);
    InlierRatio_outlier = InlierRatio_all(ind);

%% inliers
    ind = Knn_ind_GT(:, 1)>0;
    InlierNum_inlier = InlierNum_all(ind);
    InlierRatio_inlier = InlierRatio_all(ind);
    
%% all
    InlierRatio = InlierRatio_all;


InlierRatio = mean(InlierRatio);
InlierRatio_inlier = mean(InlierRatio_inlier); 
InlierRatio_outlier = mean(InlierRatio_outlier);
InlierNum_all = mean(InlierNum_all);
InlierNum_inlier = mean(InlierNum_inlier); 
InlierNum_outlier = mean(InlierNum_outlier);