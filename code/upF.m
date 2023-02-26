function  [F, eigenvalue] = upF(J_hat,K)
      
%       V is the eigenvector , D is the eigenvalue
        [V, D] = eig(J_hat);
        D = diag(D);
%       D is sorted in descending order, I is the sort index
        [D, I] = sort(D, 'descend');
        if K > length(D)
           K = length(D);
         end
         eigenvalue = D(1 : K);
            F = V(:, I(1 : K));
       
end
 
