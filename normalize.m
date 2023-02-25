function x = normalize(y)

   for i=1:size(y,2)
        x(:,i) = y(:,i)/sqrt(y(:,i)'*y(:,i));
    end
    x(find(isnan(x)==1)) = 0;
end