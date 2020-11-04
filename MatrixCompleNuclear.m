function [M,fvals,gnorms] = MatrixCompleNuclear(A,P,M,lambda,maxiter,tol)
% A : the matrix being completed
[n,d] = size(A);

f = @(M) norm(P.*(A-M))^2/2 + lambda*(norm(svd(M),1))^2/2; 
grad = @(U,S,V) P.*(A-U*S*V') + lambda * U(:,1:d)*V';
fvals = zeros(maxiter,1);
gnorms = zeros(maxiter,1);

for i = 1 : maxiter
    fvals(i) = f(M);
    [U,S,V] = svd(M + P.*(A-M));
    Snew = Shrink(S, lambda);
    M = U*Snew*V';
    gnorms(i) = norm(grad(U,Snew,V));
    if gnorms(i) < tol
        gnorms(i+1:end) = [];
        fvals(i+1:end) = [];
        fprintf('solved with %d iterations\n',i);
        break
    end
end

end


function S = Shrink(S, lambda)
    % S: diagnal matrix (sigular values)
    C = S>=lambda;
    S(C) = S(C) - lambda;
end