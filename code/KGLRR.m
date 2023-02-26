clear all
clc
%% Real Data result
%单细胞数据集路径
addpath(genpath(pwd))
dataset = {'Treutlin','Deng','Ting','Engel','Grover','Darmanis','Pollen','Goolam','Kolod'};

for num =1:1
    load(['Data_' dataset{num}]);
    n_space = length(unique(true_labs));
    K = n_space;
    [X,] = FilterGenesZero(in_X);

    %L2-normal
    X = normalize(X');
    result= zeros(2,9);
    NMI= zeros(1,30);
    ARI= zeros(1,30);
    j =1;

    la = 2^(-2.2);
    [ F_hat,~,err]= lrr_relaxed(X,K,la);
    for iiii= 1:100
        [grps] = kmeans(F_hat, n_space);
        NMI(iiii)=Cal_NMI(true_labs, grps);
        ARI(iiii)=Cal_ARI(true_labs, grps);

    end
    %NMI
    MAX_NMI=max(NMI);% max_NMI
    MIN_NMI=min(NMI);% min_NMI
    MEAN_NMI=sum(NMI)/100; %average_NMI
    VAR_NMI=var(NMI);

    %ARI
    MAX_ARI=max(ARI);% max_NMI
    MIN_ARI=min(ARI);% min_NMI
    MEAN_ARI=sum(ARI)/100; %average_NMI
    VAR_ARI=var(ARI);%variance

    result(1,j) = MEAN_NMI;
    result(2,j)=MEAN_ARI;
    result(3,j) = MAX_NMI;
    result(4,j)= MAX_ARI;
    result(5,j)=MIN_NMI;
    result(6,j) = MIN_ARI ;
    result(7,j)=VAR_NMI;
    result(8,j) = VAR_ARI;
end