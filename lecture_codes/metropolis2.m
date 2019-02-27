% basic Metropolis algorithm in 2D

clear all
close all

% sufficient statistics of the data
ybar = 0;
s = 1;
N = 10;

% number of MCMC steps
n_mc = 1e4;

% objective function is posterior of (mu, sigma^2) for y_i ~ N(mu, sigma^2)
% i = 1...N
% ybar = sample mean (y)
% s = sample standard deviation (y)

% pars(1) = mu
% pars(2) = sigma^2

posterior = @(pars) (pars(2) > 0) * pars(2)^-(N+2) * exp(-1*((N-1)*s^2 + N*(ybar-pars(1))^2)/(2*pars(2)^2));

% determine proposal scales in the two parameter dimensions
tau_mu = 0.6;
tau_sigma2 = 0.1;

% construct a covariance matrix for the proposals
jump_cov = [tau_mu^2, 0; 0, tau_sigma2^2];

jump_cov = jump_cov / 1;

% compute lower cholesky demposition
lchol = chol(jump_cov,'lower');

pars_mc = zeros(n_mc,2);

pars_mc(1,:) = [5.67,2];

acceptances = 0;

for i=1:n_mc
    
    pars_curr = pars_mc(i,:)';
   
    % generate a proposal from a multivariate normal distribution
    % using built in multivariate Gaussian function
    %pars_prop =  mvnrnd(pars_curr',jump_cov)';
    
    % generate a proposal from a multivariate normal distribution
    % using precomputed Cholesky decomposition
    pars_prop = pars_curr + lchol*randn(2,1);
    
    post_curr = posterior(pars_curr);
    
    post_prop = posterior(pars_prop);

    % calculate Metropolis ration
    r = post_prop / post_curr;
    
    % pick a uniform random variate between 0 and 1
    u = rand;
    
    % accept proposal with probability min(1,r)
    if u < r
        pars_mc(i+1,:) = pars_prop;
        acceptances = acceptances + 1;
    else
        pars_mc(i+1,:) = pars_curr;
    end
    
end

acc_ratio = acceptances/n_mc

figure(1)
subplot(2,1,1)
plot(pars_mc(:,1),'LineWidth',2)
ylabel('\mu')
title(['Markov Chain Trace Plot: N_{mc} = ' num2str(n_mc,'%.0f')])

subplot(2,1,2)
plot(pars_mc(:,2),'LineWidth',2)
xlabel('Chain step')
ylabel('\sigma^2')

%%

figure(2)
set(gca,'FontSize',16)
plot(pars_mc(:,1),pars_mc(:,2))
ylabel('\sigma^2','FontSize',16)
xlabel('\mu','FontSize',16)
title(['MCMC 2D Trace Plot: acc ratio = ' num2str(acc_ratio,'%.3f')])

pars_mc = pars_mc(n_mc/2 : end,:);

%% plot histogram

mu_grid = (-2:0.01:2)';
s2_grid = (0.5:0.01:3.25)';

figure(4)
subplot(1,2,1)
h1=histogram(pars_mc(:,1),'Normalization','pdf');
xlabel('\mu','FontSize',18)
ylabel('Posterior P(\mu | y)','FontSize',18)
hold on

% theoretical marginal posterior pdf P(mu | y) is a t distribution
plot(mu_grid,tpdf((ybar-mu_grid)*sqrt(N)/s,N-1)*sqrt(N)/s,'LineWidth',2)
hold off

subplot(1,2,2)
h2=histogram(pars_mc(:,2),'Normalization','pdf');
xlabel('\sigma^2','FontSize',18)
ylabel('Posterior P(\sigma | y)','FontSize',18)
hold on
% theoretical marginal posterior pdf P(\sigma^2 | y) 
% is a scaled inv chi^2 distribution
margpdf_s2 = (s2_grid).^-(N+1)/2 .* exp(-(N-1)*s^2 ./ (2*s2_grid.^2));
% normalise
margpdf_s2 = margpdf_s2 /(0.01*sum(margpdf_s2));
plot(s2_grid,margpdf_s2,'LineWidth',2)

hold off

%%

figure(3)
plot(pars_mc(:,1), pars_mc(:,2),'.','Color',[0.5,0.5,0.5])
xlabel('\mu','FontSize',18)
ylabel('\sigma^2','FontSize',18)
title('Posterior P(\mu, \sigma | y)','FontSize',18)
hold on

[X,Y] = meshgrid(mu_grid,s2_grid);
%post_grid = Y.^(-(N+2)/2) .* exp(-0.5*( (N-1)*s^2 + N*(ybar-X).^2)./Y);
logpost_grid = -0.5*(N+2) * Y - 0.5*( (N-1)*s^2 + N*(ybar-X).^2)./Y;
logpost_grid = logpost_grid - max(max(logpost_grid));

% evaluate the theoretical posterior contours
contour(X,Y,exp(logpost_grid),[0.0001,0.001,0.01,0.05:0.1:0.95],'LineWidth',2)

hold off



