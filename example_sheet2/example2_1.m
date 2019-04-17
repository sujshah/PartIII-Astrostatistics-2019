clear all;
close all;

%font size
fs = 14

%file = 'example_sheet2_prb1_data.txt';
file = 'ex2p1_data.txt'

fid = fopen(file);

table = textscan(fid,'%f %f %f %f','CommentStyle','#');

xs = table{1};
xerrs = table{2};
ys = table{3};
yerrs = table{4};

N = length(xs);

figure(1)
errorbar_xy2(xs,ys,xerrs,yerrs,'MarkerSize',20)

set(gca,'FontSize',fs)

xlabel('Measured x','FontSize',fs)
ylabel('Measured y','FontSize',fs)
title(['Data : N = ' num2str(N)])

xlim([-4,2])
ylim([-2,6])

%%

% run likelihood maximization

[params,ee,L,fisher] = mle_regress_xy_gx1(ys,yerrs,xs,xerrs)

alpha_fit = params(1)
beta_fit = params(2)
sigma_fit = params(3)

alpha_err = ee(1)
beta_err = ee(2)
sigma_err = ee(3)

xis_range = (-4:0.1:2)';

true_alpha = 3;
true_beta = 1;

figure(2)

errorbar_xy2(xs,ys,xerrs,yerrs,'MarkerSize',20);
hold on
set(gca,'FontSize',fs)



xlabel('Measured x','FontSize',fs)
ylabel('Measured y','FontSize',fs)
title(['Data : N = ' num2str(N)])

h1 = plot(xis_range,true_alpha + true_beta*xis_range,'-k','LineWidth',2);
hfit = plot(xis_range,alpha_fit + beta_fit*xis_range,'--k','LineWidth',3);
hold on
title('MLE fit')
legend([h1,hfit],{['True \alpha = ' num2str(true_alpha,'%.2f') ', \beta = ' num2str(true_beta,'%.2f')],...
    ['MLE \alpha = ' num2str(alpha_fit,'%.2f'), '\pm ' num2str(ee(1),'%.2f') ', \beta = ' num2str(beta_fit,'%.2f') '\pm ' num2str(ee(2),'%.2f')]'},'Location','NorthWest','FontSize',fs,'Box','Off');
hold off

%% run likelihood maximization with \tau --> infty

% in practice could use same likelihood code but fixing
% tau = 1e6 very large

[params,ee,L] = mle_regress_xy_gx1(ys,yerrs,xs,xerrs,[nan,nan,nan,nan,1e3]);

alpha_fit = params(1)
beta_fit = params(2)
sigma_fit = params(3)

alpha_err = ee(1)
beta_err = ee(2)
sigma_err = ee(3)

xis_range = (-4:0.1:2)';

true_alpha = 3;
true_beta = 1;

figure(3)

errorbar_xy2(xs,ys,xerrs,yerrs,'MarkerSize',20)

set(gca,'FontSize',fs)

xlabel('Measured x','FontSize',fs)
ylabel('Measured y','FontSize',fs)
title(['Data : N = ' num2str(N)])
hold on

h2 = plot(xis_range,true_alpha + true_beta*xis_range,'-k','LineWidth',2);
hfit2 = plot(xis_range,alpha_fit + beta_fit*xis_range,'--k','LineWidth',3);
title('MLE with \tau -> \infty fit')
hold on
legend([h2,hfit2],{['True \alpha = ' num2str(true_alpha,'%.2f') ', \beta = ' num2str(true_beta,'%.2f')],...
    ['MLE Flat \alpha = ' num2str(alpha_fit,'%.2f'), '\pm ' num2str(ee(1),'%.2f') ', \beta = ' num2str(beta_fit,'%.2f') '\pm ' num2str(ee(2),'%.3f')]'},'Location','NorthWest','FontSize',fs,'Box','Off');

%% RSS

rss = @(theta) sum( (ys - theta(1) - theta(2)*xs).^2);

start = [4;2];

[out,fval,exitflag,output,lambda,grad,hessian] = fmincon(rss,start,[],[],[],[],[],[],[]);

alpha_fit = out(1)
beta_fit = out(2)

X =  [ones(N,1), xs];
cov_err = inv(X'*X)*fval/(N-2);
ee = sqrt(diag(cov_err))

figure(4)

errorbar_xy2(xs,ys,xerrs,yerrs,'MarkerSize',20)

set(gca,'FontSize',fs)

xlabel('Measured x','FontSize',fs)
ylabel('Measured y','FontSize',fs)
title(['Data : N = ' num2str(N)])
hold on

h3 = plot(xis_range,true_alpha + true_beta*xis_range,'-k','LineWidth',2);
hfit3 = plot(xis_range,alpha_fit + beta_fit*xis_range,'--k','LineWidth',3);
title('min RSS')
hold on
legend([h3,hfit3],{['True \alpha = ' num2str(true_alpha,'%.2f') ', \beta = ' num2str(true_beta,'%.2f')],...
    ['OLS \alpha = ' num2str(alpha_fit,'%.2f'), '\pm ' num2str(ee(1),'%.2f') ', \beta = ' num2str(beta_fit,'%.2f') '\pm ' num2str(ee(2),'%.3f')]'},'Location','NorthWest','FontSize',fs,'Box','Off');

%% min chi^2

obj = @(theta) sum( (ys - theta(1) - theta(2)*xs).^2 ./yerrs.^2 );

start = [4;2];

[out,fval,exitflag,output,lambda,grad,hessian] = fmincon(obj,start,[],[],[],[],[],[],[]);

min_chi2 = fval
alpha_fit = out(1)
beta_fit = out(2)

X =  [ones(N,1), xs];
Wmat = diag(yerrs.^2);
Varb = inv(X' * inv(Wmat) * X)
ee = sqrt(diag(cov_err))

figure(5)

errorbar_xy2(xs,ys,xerrs,yerrs,'MarkerSize',20)

set(gca,'FontSize',fs)

xlabel('Measured x','FontSize',fs)
ylabel('Measured y','FontSize',fs)
title(['Data : N = ' num2str(N)])
hold on

h3 = plot(xis_range,true_alpha + true_beta*xis_range,'-k','LineWidth',2);
hfit3 = plot(xis_range,alpha_fit + beta_fit*xis_range,'--k','LineWidth',3);
title('min \chi^2')
hold on
legend([h3,hfit3],{['True \alpha = ' num2str(true_alpha,'%.2f') ', \beta = ' num2str(true_beta,'%.2f')],...
    ['WLS \alpha = ' num2str(alpha_fit,'%.2f'), '\pm ' num2str(ee(1),'%.2f') ', \beta = ' num2str(beta_fit,'%.2f') '\pm ' num2str(ee(2),'%.3f')]'},'Location','NorthWest','FontSize',fs,'Box','Off');

%% min EXY

obj = @(theta) sum( (ys - theta(1) - theta(2)*xs).^2 ./(yerrs.^2 +theta(2)^2 * xerrs.^2) );

start = [4;2];

[out,fval,exitflag,output,lambda,grad,hessian] = fmincon(obj,start,[],[],[],[],[],[],[]);

min_EXY = fval
alpha_fit = out(1)
beta_fit = out(2)

figure(6)

errorbar_xy2(xs,ys,xerrs,yerrs,'MarkerSize',fs)

set(gca,'FontSize',18)

xlabel('Measured x','FontSize',fs)
ylabel('Measured y','FontSize',fs)
title(['Data : N = ' num2str(N)])
hold on

h3 = plot(xis_range,true_alpha + true_beta*xis_range,'-k','LineWidth',2);
hfit3 = plot(xis_range,alpha_fit + beta_fit*xis_range,'--k','LineWidth',3);
title('min EXY')
hold on
legend([h3,hfit3],{['True \alpha = ' num2str(true_alpha,'%.2f') ', \beta = ' num2str(true_beta,'%.2f')],...
    ['EXY \alpha = ' num2str(alpha_fit,'%.2f') ', \beta = ' num2str(beta_fit,'%.2f') ]'},'Location','NorthWest','FontSize',fs,'Box','Off');

