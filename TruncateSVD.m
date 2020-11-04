function [Ak] = TruncateSVD(A,k)
    [U,S,V] = svd(A);
    Ak = U(:,1:k)*S(1:k,1:k)*V(:,1:k)';
end