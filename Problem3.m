% Problem 3
close all
[M,y] = readdata();
runs = 100;

[MU,MS,MV] = svd(M);

rel_err = zeros(8,9);
abs_err = zeros(8,9);

for a = 1 : 8
    for k = 2 : 10
        c = a*k;
        r = a*k;
        
        Mk = MU(:,1:k)*MS(1:k,1:k)*MV(:,1:k)';
        errors = zeros(runs,1);
        abs_errors = 0;
        for run = 1 : runs
            [C,U,R] = CUR(M,c,r,k,MU,MV);
            errors(run) = norm(M-C*U*R,'fro');
        end
        rel_err(a,k-1) = mean(errors)/norm(M-Mk,'fro');
        abs_err(a,k-1) = mean(errors);
        fprintf('a=%d, rank %d, mean error %d / %d \n',a, k, abs_err(a,k-1), rel_err(a,k-1))
    end
    
end

fsz = 16;clf;
figure(1);
hold on;
grid;
for a = 1 : 8
    plot((2:10)',rel_err(a,:),'Linewidth',2,'DisplayName',['a',num2str(a)]);
end
set(gca,'Fontsize',fsz);
legend('show');
xlabel('k','Fontsize',fsz);
ylabel('relative error','Fontsize',fsz);
title('Relative Error')
filename = 'figs/problem3-1.png';
saveas(gcf,filename)

figure(2);
hold on;
grid;
for a = 1 : 8
    plot((2:10)',abs_err(a,:),'Linewidth',2,'DisplayName',['a',num2str(a)]);
end
set(gca,'Fontsize',fsz);
legend('show');
xlabel('k','Fontsize',fsz);
ylabel('absolute error','Fontsize',fsz);
title('Absolute Error')
filename = 'figs/problem3-2.png';
saveas(gcf,filename)

