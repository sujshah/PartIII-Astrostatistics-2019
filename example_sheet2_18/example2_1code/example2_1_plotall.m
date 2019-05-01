clear all;
close all;

file = 'example_sheet2_prb1_data.txt';

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

xis_range = (-4:0.1:1)';

true_alpha = 3;
true_beta = 1;

alpha_mle = 3.26;
beta_mle  = 0.90;

alpha_flat = 2.75;
beta_flat  = 0.46;

alpha_ols = 2.75;
beta_ols  = 0.46;

alpha_exy = 4.06;
beta_exy  = 1.58;

alpha_post = 3.45;
beta_post  = 1.06;

htrue = plot(xis_range,true_alpha + true_beta*xis_range,'-k','LineWidth',2);
hmle = plot(xis_range,alpha_mle + beta_mle*xis_range,'--k','LineWidth',3);
hols = plot(xis_range,alpha_ols + beta_ols*xis_range,'--b','LineWidth',2);
hexy = plot(xis_range,alpha_exy + beta_exy*xis_range,'--r','LineWidth',2);
hpost = plot(xis_range,alpha_post + beta_post*xis_range,'--m','LineWidth',3);

legend([htrue,hmle,hols,hexy,hpost],{'True \beta = 1','MLE \beta = 0.90 \pm 1.2','OLS/\chi^2/Flat \beta = 0.46','EXY \beta = 1.58','Posterior \beta = 1.06 \pm 0.47'},'Location','NorthWest','Box','Off')

xlim([-4,1])
ylim([0,4])

