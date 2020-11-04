function [L,R,fvals] = KMeans(A,k,maxiter)
% A : the matrix being decomposed
[n,d] = size(A);

R = A(1:k,:); % initialize the mean of clusters as the first k rows
L = zeros(n,k);
L(1:k,:) = eye(k);

f = @(L,R) norm(A-L*R, 'fro')^2/2;
fvals = zeros(maxiter,1);

for i = 1 : maxiter
    flag = 0;
    fvals(i) = f(L,R);
    
    for j = 1 : n
        ind = FindNearest(A(j,:),R,k);
        if find(L(j,ind)==1)
%             fprintf("not change\n")
            continue
        else
            cur_cluster = L(j,:)==1;
            L(j,cur_cluster) = 0;
            L(j,ind(1)) = 1;
            flag = 1;
        end
    end
    
    
    if flag == 0
        fprintf('solved with %d iterations\n',i);
        break
    end
    
    % recalculate means
    for t = 1 : k
        in_cluster = L(:,t)==1;
        if max(in_cluster) == 0
            randi = randperm(n);
            R(t,:) = A(randi(1),:);
        else
            new_mean = mean(A(in_cluster,:));
            R(t,:) = new_mean;
        end
    end
end

end

function ind = FindNearest(row,clusters,k)
    distances = zeros(k,1);
    for i = 1 : k
        distances(i) = norm(row' - clusters(i,:)');
    end
    ind = find(distances==min(distances));
end