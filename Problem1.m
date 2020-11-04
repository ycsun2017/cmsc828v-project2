% Problem 1: NMF
A = csvread('MovieRankings36.csv');
[n,d] = size(A);
P = A~=0;
%% first complete the matrix
lambda = 4;
k = 5;
X = rand(n,k);
Y = rand(d,k);
[X,Y,~,~] = MatrixCompleAlter(A,P,X,Y,lambda,1000,1e-2,1e-5);
fprintf('lambda %d, final error %d\n',lambda, norm(P.*(A-X*Y')));
M = X*Y';

%% K-means algorithm
k = 5;
[L,R,fs] = KMeans(M,k,1000);
fprintf('k=%d means error %d\n', k, norm(A-L*R));

%% Projected gradient descent
k = 5;
W = rand(n,k);
H = rand(k,d);
[pgd_W,pgd_H,pgd_fs,pgd_gs] = PGD(M,W,H,1000,1e-3,1e-3);
fprintf('PGD: rank %d, final error %d\n',k,norm(M-pgd_W*pgd_H));

%% Lee-Seung scheme
[ls_W,ls_H,ls_fs,ls_gs] = LeeSeung(M,W,H,1000,1e-3);
fprintf('Lee-Seung: rank %d, final error %d\n',k,norm(M-ls_W*ls_H));

%% PGD + Lee-Seung 
[pls_W,pls_H,pls_fs,pls_gs] = PGDLS(M,W,H,500,1000,1e-3,1e-3);
fprintf('PGD + Lee-Seung: rank %d, final error %d\n',k,norm(M-pls_W*pls_H));

%% Plotting function values and gradient norms
fsz = 16;
figure;clf;
subplot(2,1,1);
hold on;
grid;
plot((1:length(pgd_fs))',pgd_fs,'Linewidth',2,'Marker','.','Markersize',20);
plot((1:length(ls_fs))',ls_fs,'Linewidth',2,'Marker','.','Markersize',20);
plot((1:length(pls_fs))',pls_fs,'Linewidth',2,'Marker','.','Markersize',20);
legend('PGD','Lee-Seung','PGD+Lee-Seung');
set(gca,'Fontsize',fsz);
set(gca, 'YScale', 'log')
xlabel('k','Fontsize',fsz);
ylabel('f','Fontsize',fsz);
subplot(2,1,2);
hold on;
grid;
plot((1:length(pgd_gs))',pgd_gs,'Linewidth',2,'Marker','.','Markersize',20);
plot((1:length(ls_gs))',ls_gs,'Linewidth',2,'Marker','.','Markersize',20);
plot((1:length(pls_gs))',pls_gs,'Linewidth',2,'Marker','.','Markersize',20);
legend('PGD','Lee-Seung','PGD+Lee-Seung');
set(gca,'Fontsize',fsz);
set(gca, 'YScale', 'log')
xlabel('k','Fontsize',fsz);
ylabel('||g||','Fontsize',fsz);