function [W,H,fvals,gnorms] = PGDLS(A,W,H,maxiter1,maxiter2,stepsize,tol)
% Projected gradient descent + Lee-Seung
% A : the matrix being decomposed
maxiter = maxiter1 + maxiter2;
f = @(W,H) norm(A-W*H, 'fro')^2/2;
R = @(W,H) A - W*H;
gradw = @(W,H) R(W,H)*H';
gradh = @(W,H) W'*R(W,H);

fvals = zeros(maxiter,1);
gnorms = zeros(maxiter,1);

% start with PGD
for i = 1 : maxiter1
    fvals(i) = f(W,H);
    Wnew = Proj(W+stepsize*R(W,H)*H');
    Hnew = Proj(H+stepsize*W'*R(W,H));
    W = Wnew;
    H = Hnew;
    gnorms(i) = norm(gradw(W,H))+norm(gradh(W,H));
    if gnorms(i) < tol
        gnorms(i+1:end) = [];
        fvals(i+1:end) = [];
        fprintf('solved with %d iterations\n',i);
        break
    end
end

% continue with Lee-Seung
for i = maxiter1+1 : maxiter2
    fvals(i) = f(W,H);
    W = W./(W*H*H').*(A*H');
    H = H./(W'*W*H).*(W'*A);
    gnorms(i) = norm(gradw(W,H))+norm(gradh(W,H));
    if gnorms(i) < tol
        gnorms(i+1:end) = [];
        fvals(i+1:end) = [];
        fprintf('solved with %d iterations\n',i);
        break
    end
end

end

function M = Proj(M)
    % element-wise projection to nonnegative matrix set
    C = M<0;
    M(C) = 0;
end