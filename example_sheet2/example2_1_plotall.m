clear all;
close all;

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
hold on
set(gca,'FontSize',18)

xlabel('Measured x','FontSize',18)
ylabel('Measured y','FontSize',18)
title(['Data : N = ' num2str(N)])

xis_range = (-4:0.1:2)';

true_alpha = 3;
true_beta = 1;

alpha_mle = 3.05;
beta_mle  = 1.1;

alpha_flat = 2.63;
beta_flat  = 0.57;

alpha_ols = 2.63;
beta_ols  = 0.55;

alpha_wls = 2.52;
beta_wls = 0.51;

alpha_exy = 3.13;
beta_exy  = 1.33;

alpha_post = 3.04;
beta_post  = 1.087;

htrue = plot(xis_range,true_alpha + true_beta*xis_range,'-k','LineWidth',2);
hmle = plot(xis_range,alpha_mle + beta_mle*xis_range,'--k','LineWidth',3);
hols = plot(xis_range,alpha_ols + beta_ols*xis_range,'--b','LineWidth',2);
hwls = plot(xis_range,alpha_wls + beta_wls*xis_range,'--b','LineWidth',2);
hexy = plot(xis_range,alpha_exy + beta_exy*xis_range,'--r','LineWidth',2);
hpost = plot(xis_range,alpha_post + beta_post*xis_range,'--m','LineWidth',3);

legend([htrue,hmle,hols,hwls,hexy,hpost],{'True','MLE','OLS/MLE Flat','min chi^2','EXY','Posterior'},'Location','NorthWest','Box','Off')

xlim([-4,2])
ylim([-2,6])

