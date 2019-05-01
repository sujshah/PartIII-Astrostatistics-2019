clear all
close all;

n = 1000;

% set the true values and variances
Mbar_true = -18;
sigma_int_true = 0.1;

sigma_ms = 0.05*ones(n,1);
sigma_mus = 0.1*ones(n,1);
sigmas = sqrt(sigma_ms.^2 + sigma_mus.^2);

% generate data;

Ms_true = Mbar_true + randn(n,1)*sigma_int_true;

Mobs = Ms_true + randn(n,1).*sigmas;

% plot the data
figure(1);

edges = -18.7:0.1:-17;

htrue=histogram(Ms_true,edges);
hold on
hobs=histogram(Mobs,edges);
hold on;

hold off
set(gca,'FontSize',20)

xlabel('Absolute Magnitude','FontSize',20)
ylabel('Number','FontSize',20)
legend([htrue,hobs],{'Latent (True) Abs Mags','Observed Abs Mag'},'FontSize',18,'Location','NorthEast')
title(['Sample STD(true) = ' num2str(std(Ms_true),'%.3f') ', Sample STD(obs) = ' num2str(std(Mobs),'%.3f')],'FontSize',18)

xlim([min(edges),max(edges)])
ymax = 15+max(htrue.Values);
ylim([0,ymax])

text(-18.6,ymax-50,['Mbar_{true} = ' num2str(Mbar_true,'%.2f')],'FontSize',20)
text(-18.6,ymax-75,['\sigma_{int,true} = ' num2str(sigma_int_true,'%.2f')],'FontSize',20)
text(-18.6,ymax-25,['N = ' num2str(n)],'FontSize',20)

%% fit data

disp('Finding Maximum Likelihood Solution')

% a quick function of theta = (mu, sigma_int) to evaluate the negative log likelihood 
neglogL = @(theta) 0.5* sum( (Mobs-theta(1)).^2 ./(sigmas.^2 + theta(2)^2) + log(2*pi*(sigmas.^2 + theta(2)^2)) );

% set some options
options=optimset('Display','iter','UseParallel',0,'Algorithm','active-set');
   
init = [mean(Mobs),std(Mobs)];

lb = [-Inf;0];
ub = [Inf;Inf];

[bestfit,fval,exitflag,output,lambda,grad,hessian]= fmincon(neglogL,init,[],[],[],[],lb,ub,[],options);

mu_fit = bestfit(1)
sigma_fit = bestfit(2)

cov_unc = inv(hessian)

param_errs = sqrt(diag(cov_unc));

mu_err = param_errs(1)
sigma_err = param_errs(2)

disp('-----Maximum Likelihood Results-----')
disp(['mu_hat = ' num2str(mu_fit,'%.3f') ' +/- ' num2str(mu_err,'%.3f') ' mag'])
disp(['sigma_int_hat = ' num2str(sigma_fit,'%.3f') ' +/- ' num2str(sigma_err,'%.3f') ' mag'])
disp('------------------------------------')

%% evaluate likelihood contours

mu_grid = (mu_fit-3.5*mu_err):(mu_err/10):(mu_fit+3.5*mu_err);
sigma_grid = (sigma_fit-3.5*mu_err):(mu_err/10):(sigma_fit+3.5*mu_err);

[X,Y] = meshgrid(mu_grid, sigma_grid);
loglkhd_grid = X*0;

logL = @(theta) -0.5* sum( (Mobs-theta(1)).^2 ./(sigmas.^2 + theta(2)^2) + log(2*pi*(sigmas.^2 + theta(2)^2)) );

disp('Evaluating Likelihood Grid...');
for i=1:length(mu_grid)
    for j=1:length(sigma_grid)
        
        loglkhd_grid(i,j) = logL([mu_grid(i),sigma_grid(j)]);
    end
end
disp('Done!');

% this just normalizes the likelihood, so the value at the peak is 1.
lkhd_grid = exp(loglkhd_grid - max(max(loglkhd_grid)));

figure(2)
set(gca,'FontSize',18)
contour(X,Y,lkhd_grid,[0.05:0.1:0.95])
hold on

% rule of thumb, when the likelihood is Gaussian in the paramters enough
% the 2D 68% and 95% contours can be found by finding the contour levels 
% that are 0.32 and 0.05 times the max lkhd value.s
% see Gelman BDA, Sec 4.1 (pg 85)

h=contour(X,Y,lkhd_grid,[0.05, 0.32],'-k','LineWidth',2);
plot(mu_fit,sigma_fit,'xk','MarkerSize',25)

xlabel('\mu','FontSize',20)
ylabel('\sigma','FontSize',20)
title('Likelihood Contours and 95%, 68% Confidence Regions','FontSize',18)

xlim([mu_fit - 3.25*0.015,mu_fit + 3.25*0.015])
ylim([sigma_fit - 3.25*0.015,sigma_fit + 3.25*0.015])
