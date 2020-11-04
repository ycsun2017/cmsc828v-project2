function [C,U,R] = CUR(M,c,r,k,MU,MV)
% CUR decomposition
    C = ColumnSelect(M,k,c,MV);
    R = ColumnSelect(M',k,r,MU)';
    U = pinv(C)*M*pinv(R);
end

 
function [CorR] = ColumnSelect(A,k,c,V)
    Vk = V(:,1:k);
    p = c * sum(Vk.^2,2) / k; 
    p(p > 1) = 1;
    CorR = A(:,(rand(size(p))<=p));
end