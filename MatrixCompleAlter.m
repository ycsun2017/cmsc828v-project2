function [X,Y,fvals,gnorms] = MatrixCompleAlter(A,P,X,Y,lambda,maxiter,stepsize,tol)
% A : the matrix being completed
% [n,d] = size(A);
J = 3;
f = @(X,Y) norm(P.*(A-X*Y'))^2/2 + lambda*(norm(X)^2+norm(Y)^2)/2; 
gradx = @(X,Y) - P.*(A-X*Y')*Y + lambda*X;
grady = @(X,Y) - P'.*(A-X*Y')'*X + lambda*Y;

fvals = zeros(maxiter,1);
gnorms = zeros(maxiter,1);

for i = 1 : maxiter
    fvals(i) = f(X,Y);
    % update X
    for j = 1 : J
        X = X - stepsize * gradx(X,Y);
    end
    % update Y
    for j = 1 : J
        Y = Y - stepsize * grady(X,Y);
    end
    gnorms(i) = norm(gradx(X,Y)) + norm(gradx(X,Y));
    if gnorms(i) < tol
        gnorms(i+1:end) = [];
        fvals(i+1:end) = [];
        fprintf('solved with %d iterations\n',i);
        break
    end
end

end