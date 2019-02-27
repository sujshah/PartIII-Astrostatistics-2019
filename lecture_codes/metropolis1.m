% basic Metropolis algorithm

clear all
close all

ybar = 0
N = 10
sigma = 1

posterior = @(mu) normpdf(mu,ybar,sigma/sqrt(N));

n_mc = 2000;
tau = 0.6;

mus_mc = zeros(n_mc,1);

mus_mc(1) = 5.67;

acceptances = 0;

for i=1:n_mc
    
    mu_curr = mus_mc(i);
    
    mu_prop = mu_curr + randn*tau;
    
    post_curr = posterior(mu_curr);
    
    post_prop = posterior(mu_prop);
    
    r = post_prop / post_curr;
    
    % pick a uniform random variate between 0 and 1
    u = rand;
    
    if rand < r
        mus_mc(i+1) = mu_prop;
        acceptances = acceptances + 1;
    else
        mus_mc(i+1) = mu_curr;
    end
    
end

acc_ratio = acceptances/n_mc

%%

figure(1)
subplot1(2,1,'Gap',[0.03,0.03])
subplot1(1)
plot(mus_mc,'LineWidth',2)
ylabel('Parameter Value \mu','FontSize',20)
%xlabel('Chain step','FontSize',20)
title(['Markov Chain Trace Plot: acc ratio = ' num2str(acc_ratio,'%.2f')],'FontSize',18)
set(gca,'FontSize',18)
subplot1(2)
plot(log(posterior(mus_mc)),'LineWidth',2)
ylabel('Log Posterior Log P(\mu | y)','FontSize',18)
xlabel('Chain step','FontSize',20)
set(gca,'FontSize',18)
%% plot histogram

mus_mc = mus_mc(1 : end);
mus_mc = mus_mc(n_mc/2 : 2:  end);

figure(2)

h=histogram(mus_mc,-2.125:0.25:2.125,'Normalization','pdf');
hold on
x = -2:0.01:2;
h2=plot(x, normpdf(x,0,1/sqrt(N)),'LineWidth',2);
xlabel('\mu','FontSize',18)
ylabel('Posterior P(\mu | y)','FontSize',18)
legend({'Posterior Histogram','Analytic Posterior'})
title(['N_{mc} = ' num2str(n_mc,'%.0f')],'FontSize',18)
set(gca,'FontSize',20)

