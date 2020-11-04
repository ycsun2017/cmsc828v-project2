% Problem 2: Matrix Completion
A = csvread('MovieRankings36.csv');

% parameters
k = 10;
[n,d] = size(A);
P = A~=0;
errors = zeros(10,1);
type = "nuclear";

for i = 1 : 10
    lambda = i*0.1;
    
    if type == "alter"
        % Low-rank factorization with alternating iteration
        X = rand(n,k);
        Y = rand(d,k);
        fprintf('lambda %d, initial error %d\n',lambda, norm(P.*(A-X*Y')));
        [X,Y,fs,gs] = MatrixCompleAlter(A,P,X,Y,lambda,1000,1e-2,1e-5);
        fprintf('lambda %d, final error %d\n',lambda, norm(P.*(A-X*Y')));
        errors(i) = norm(P.*(A-X*Y'));
    end
    
    if type == "nuclear"
        % Nuclear norm trick
        M = rand(n,d);
        fprintf('lambda %d, initial error %d\n',lambda, norm(P.*(A-M)));
        [M,fs,gs] = MatrixCompleNuclear(A,P,M,lambda,1000,1e-5);
        fprintf('lambda %d, final error %d\n',lambda, norm(P.*(A-M)));
        errors(i) = norm(P.*(A-M));
    end
    
end

%% Plotting function values and gradient norms
% fsz = 16;
% figure;clf;
% subplot(2,1,1);
% hold on;
% grid;
% niter = length(fs);
% plot((0:niter-1)',fs,'Linewidth',2);
% set(gca,'Fontsize',fsz);
% xlabel('k','Fontsize',fsz);
% ylabel('f','Fontsize',fsz);
% subplot(2,1,2);
% hold on;
% grid;
% niter = length(gs);
% plot((0:niter-1)',gs,'Linewidth',2);
% set(gca,'Fontsize',fsz);
% xlabel('k','Fontsize',fsz);
% ylabel('||g||','Fontsize',fsz);

%% Plotting results with different lambdas
figure;clf;
hold on;
grid;
plot((0.1:0.1:1.0)',errors,'Linewidth',2);
set(gca,'Fontsize',fsz);
xlabel('lambda','Fontsize',fsz);
ylabel('error norm','Fontsize',fsz);