function [F] = zhibiao(Rnk) 
%    µ¥Î»»¯
     [a,b] = size(Rnk);
     for i =1:b
         A (1,i)= norm(Rnk(:,i));
     end
     A = repmat(A,a,1);
     F= Rnk./A; 
     B = F'*F;
end