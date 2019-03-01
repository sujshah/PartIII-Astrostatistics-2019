clear all
close all

fid = fopen('ex1_data_for_problem3.txt');
table = textscan(fid,'%f');
yobs = table{1};
N = length(yobs);

t0 = 1
t1 = 5

%% calculate naive MLE

mle_gamma = @(yy) length(yy)/sum(log(yy/t0)) + 1;
gamma_hat = mle_gamma(yobs)
C = (gamma_hat-1)/(t0^(1-gamma_hat));

% calculate bootstrap uncertainty
[bootstat,bootsam] = bootstrp(100,mle_gamma, yobs);
mle_booterr = std(bootstat)

%% plot corrected likelihood function
figure(3)
set(gca,'FontSize',20)
grid = 1:0.001:10;

log_lkhd = (-grid).*sum(log(yobs)) + N*log( (grid-1)./(t0.^(1-grid) - t1.^(1-grid)));

mle = grid(log_lkhd == max(log_lkhd));

plot(grid,exp(log_lkhd-max(log_lkhd)),'-b','LineWidth',3)
hold on
set(gca,'LineWidth',2)
log_lkhd_wrong = (-grid).*sum(log(yobs)) + N.*log(grid-1) + N*(grid-1)*log(t0);

mle_wrong = grid(log_lkhd_wrong == max(log_lkhd_wrong));


h1=plot(grid,exp(log_lkhd_wrong-max(log_lkhd_wrong)),'-r','LineWidth',3);
hold on

xlabel('\gamma (power law exponent)','FontSize',16)
ylabel('Likelihood','FontSize',16)
title(['Observed Data N = ' num2str(N)],'FontSize',14)

legend(['Correct Likelihood (MLE = ' num2str(mle) ')'],['Naive Likelihood (MLE = ' num2str(mle_wrong) ')'],'Location','NorthEast')
hold off
set(gca,'FontSize',20)

%% calculate bootstrap of corrected MLE

yboot = bootsam*0;
mles = zeros(100,1);
mle_wrongs = zeros(100,1);

for b=1:100
    yboot(:,b) = yobs(bootsam(:,b));
    
    grid = 1:0.001:10;
    
    log_lkhd = (-grid).*sum(log(yboot(:,b))) + N*log( (grid-1)./(t0.^(1-grid) - t1.^(1-grid)));
    
    mle = grid(log_lkhd == max(log_lkhd));
    mles(b) = mle;

    log_lkhd_wrong = (-grid).*sum(log(yboot(:,b))) + N*log(grid-1) + N*(grid-1)*log(t0);
    
    mle_wrong = grid(log_lkhd_wrong == max(log_lkhd_wrong));
    mle_wrongs(b) = mle_wrong;
        
end

corr_mle_booterr = std(mles)

%% generate new data

mle_sims = zeros(100,1);

for ii=1:100
    
    % generate a random sample of N from the truncated pareto using inverse
    % CDF method
    u = rand(N,1);
    y = (t1^(1-mle)- t0^(1-mle)).*u + t0^(1-mle);
    y = y.^(1/(1-mle));
    
    grid = 1:0.001:10;

    log_lkhd = (-grid).*sum(log(y)) + N*log( (grid-1)./(t0.^(1-grid) - t1.^(1-grid)));

    mle_sims(ii) = grid(log_lkhd == max(log_lkhd));
    
end

mean_mle_sims = mean(mle_sims)
std_mle_sims = std(mle_sims)


