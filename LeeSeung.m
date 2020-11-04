function [W,H,fvals,gnorms] = LeeSeung(A,W,H,maxiter,tol)
% Lee-Seung (multiplicative update scheme)
% A : the matrix being decomposed

f = @(W,H) norm(A-W*H, 'fro')^2/2;
R = @(W,H) A - W*H;
gradw = @(W,H) R(W,H)*H';
gradh = @(W,H) W'*R(W,H);

fvals = zeros(maxiter,1);
gnorms = zeros(maxiter,1);

for i = 1 : maxiter
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
