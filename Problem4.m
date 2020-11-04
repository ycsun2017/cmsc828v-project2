% text categorization
[M,y] = readdata();
[MU,MS,MV] = svd(M);
[n,d] = size(M);
type1 = 1:71;
type2 = 72:139;

% importance score
x = pinv(M) * y;
[maxx, maxxind] = sort(abs(x), 'descend');

k = 15;
newM = M(:,maxxind(1:10000));
[U,S,V] = svd(newM);
Vk2 = V(:,1:k);
pi2 = sum(Vk2.^2,2) / k;
[maxlev2, maxlevind2] = sort(pi2, 'descend');
newM5 = newM(:,maxlevind2(1:5));

centered = newM5 - ones(n,1) * sum(newM5)/n;
[U,S,V] = svd(centered);
principled = V(:,1:2);
mapped = centered * principled;

fsz = 16;clf;
figure(1);
hold on;
grid;
scatter(mapped(type1,1),mapped(type1,2),'d', 'r');
scatter(mapped(type2,1),mapped(type2,2),'o','b');

% Vk = MV(:,1:k);
% pi = sum(Vk.^2,2) / k;
% [maxlev, maxlevind] = sort(pi, 'descend');
% 
% newM = M(:,maxlevind(1:10000));
% [U,S,V] = svd(newM);

% A = sum(M(1:71,:));
% B = sum(M(72:end,:));
% C = abs(A - B);
% [maxdiff, maxdiffind] = sort(C, 'descend');
% newM = M(:,maxdiffind(1:10000));

% dists = zeros(n,d);
% for i = type1
%     for j = type2
%        dists(i,:) = (M(i,:)-M(j,:)).^2;
%     end
% end
% 
% for i = type2
%     for j = type1
%        dists(i,:) = (M(i,:)-M(j,:)).^2;
%     end
% end
% [maxdist, maxdistind] = sort(mean(dists), 'descend');

% dists2 = zeros(n,d);
% mean1 = mean(M(type1,:));
% mean2 = mean(M(type2,:));
% dists2(type1,:) = (M(type1,:) - ones(71,1) * mean2).^2;
% dists2(type2,:) = (M(type2,:) - ones(68,1) * mean1).^2;
% [maxdist2, maxdistind2] = sort(mean(dists2), 'descend');

% kk = 100;
% x = MV(:,1:kk)*inv(MS(1:kk,1:kk))*MU(:,1:kk)'*y;



% score = zeros(1,d);
% for i = 1 : d
% %     score(maxdiffind(i)) = d-i;
% %     score(maxlevind(i)) = d-i;
% %     score(maxdistind(i)) = d-i;
%     score(maxxind(i)) = d-i;
% end
% [maxsco, maxscoind] = sort(score, 'descend');



