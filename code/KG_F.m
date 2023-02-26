%%  调参专用函数
function [ F_hat,flag ]=KG_F(X,K,lambda,beta,maxIter,mu,epsilon1 , epsilon2  )
% min 1/2|| X - XS ||_F^2 + alf || J ||_* +bate tr(ZLZ') -tr(F'Z'SF)  ,beta
% addpath PROPACK;
addpath PROPACK;
n = size(X, 2);

if nargin < 7
    epsilon1 = 10^(-6);
end
if nargin < 8
    epsilon2 = 10^(-3);
end


fea = X;
options = [];
options.maxIter = 100;
options.alpha = 100;
options.error = 1e-5;
W = constructW(fea',options);

Y = X;
[mFea,nSmp]=size(Y);
DCol = full(sum(W,2));
D = spdiags(DCol,0,nSmp,nSmp);
D = full(D);
L = D - W;


S_hat = zeros(n);
Z_hat = zeros(n);

F_hat = zeros(n,K);
F_k =zeros(n,K);
Y1_hat = zeros(n);
Y2_hat = zeros(K,K);

converged = false;
flag= 0;
iter = 0;
% mu = 1e-6; % this one can be tuned

var_X = Y'*Y;
I = speye(n);
AA = inv(2* beta * (L'+ L )  + mu*I);
A = inv(var_X + mu*I);
sv = 10;
total_svd = 0;
stopCriterion = 1;
mu_max = 1e6; % this one can be tuned
rho = 1.2 ; % this one can be tuned


while ~converged
    iter = iter + 1;
    
    %         更新核范数Z
    Z_k = Z_hat;
    %          Q = (2* beta * (L'+ L )  + mu*I)*Z_hat -(mu*S_hat +S_hat *F_hat*F_hat'-Y1_hat) ;
    p = 2* beta * norm(L,2)  + mu*1;
    Q = AA*(mu*S_hat +S_hat *F_hat*F_hat'-Y1_hat) ;
    pp = (Q+Q')/2;
    
    if any(any(isnan(pp))) || any(any(isinf(pp)))
        
        disp(['第' num2str(iter) '次，Z_hat存在ANY 或者 INF值']);
        F_hat= F_k;
        flag = 1;
        break;
    end
    [Z_hat, ~] = solve_nn(pp, lambda / mu );
    %          Z_hat(Z_hat<0)=0;
    Z_hat=Z_hat-diag(diag(Z_hat));
    
    
    %       更新F
    F_k =F_hat;
    [f,d,flag]= upF(Z_hat,Z_k,K,iter);
    
    
    F_hat = zhibiao(f);
    F_hat(F_hat<0)=0;
    
    
    
    %         更新S
    S_k = S_hat  ;
    S_hat = A * (var_X + Z_hat * F_hat * F_hat' + mu * Z_hat + Y1_hat);
    S_hat(S_hat < 0) = 0 ;
    
    
    total_svd = total_svd + 1 ;
    
    %         更新Y1  Y2
    Y1_hat = Y1_hat + mu * (Z_hat - S_hat);
    Y2_hat = Y2_hat +mu*(F_hat'*F_hat - speye(K));
    
    mu = min(mu_max , rho * mu);
    
    b = norm( S_hat - S_k ,2);
    c = norm( Z_hat - Z_k ,2);
    d = norm( F_hat - F_k ,2);
    
    a = [p*c ,mu*b,mu*d ];
    a = max(a);
    
    if a <= epsilon2
        rho_hat = rho;
    else rho_hat = 1;
    end
    
    stopCriterion1 = norm( X - X * S_hat )/norm(X , 2);
    e = [p * norm(Z_hat - Z_k ,2),mu * norm(F_hat - F_k ,2), mu * norm(S_hat - S_k ,2)];
    stopCriterion2 = max(e);
    if  stopCriterion1<epsilon1 && stopCriterion2<epsilon2
        converged = 1;
    end
    
    if ~converged && iter >=maxIter
%         disp('Maximum iterations reached') ;
        converged = 1 ;
    end
    
end


