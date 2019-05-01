clear all
close all;

filename = 'example_sheet1_prb1_data.txt';

disp(['Reading in data from ' filename])

fid = fopen(filename);

table = textscan(fid,'%f %f %f %f %f','CommentStyle','#');

% unpack the data
N = length(table{1});
mhats = table{2};
sigma_ms = table{3};

muhats = table{4};
sigma_mus = table{5};

Mobs = mhats - muhats;

sigmas = sqrt(sigma_mus.^2 + sigma_ms.^2);

% plot the data
figure(1);

edges = -18.7:0.1:-17;

hobs=histogram(Mobs,edges);
hold on;

hold off
set(gca,'FontSize',20)

xlabel('Absolute Magnitude','FontSize',20)
ylabel('Number','FontSize',20)

legend([hobs],{'Observed Abs Mag'},'FontSize',18,'Location','NorthEast')
title(['N = ' num2str(N,'%.0f') ', Sample STD(obs) = ' num2str(std(Mobs),'%.3f')],'FontSize',18)

xlim([min(edges),max(edges)])
ymax = 1+max(hobs.Values);
ylim([0,ymax])

%% fit data

disp('Finding Maximum Likelihood Solution')

% a quick function of theta = (mu, sigma_int) to evaluate the negative log likelihood 
neglogL = @(theta) 0.5* sum( (Mobs-theta(1)).^2 ./(sigmas.^2 + theta(2)^2) + log(2*pi*(sigmas.^2 + theta(2)^2)) );

% set some options
options=optimset('Display','iter','UseParallel',0,'Algorithm','active-set');
   
init = [mean(Mobs),std(Mobs)];

lb = [-Inf;0];
ub = [Inf;Inf];

% optimise the function to get the MLE
% also outputs the hessian of the negative log lkhd at the MLE 
% this matrix is the observed fisher information
[bestfit,fval,exitflag,output,lambda,grad,hessian]= fmincon(neglogL,init,[],[],[],[],lb,ub,[],options);

Mbar_fit = bestfit(1)
sigma_fit = bestfit(2)

% invert the fisher matrix to get an asymptotic approximation to the
% covariance matrix
cov_unc = inv(hessian)

% the diagonals are lower bounds on the variance of the MLEs
param_errs = sqrt(diag(cov_unc));

Mbar_err = param_errs(1)
sigma_err = param_errs(2)

disp(' ')
disp('-----Maximum Likelihood Results-----')
disp(['Mbar_hat = ' num2str(Mbar_fit,'%.3f') ' +/- ' num2str(Mbar_err,'%.3f') ' mag'])
disp(['sigma_int_hat = ' num2str(sigma_fit,'%.3f') ' +/- ' num2str(sigma_err,'%.3f') ' mag'])
disp('------------------------------------')
disp(' ')

%% evaluate likelihood contours

mu_grid = (Mbar_fit-3.5*Mbar_err):(Mbar_err/10):(Mbar_fit+3.5*Mbar_err);
sigma_grid = (sigma_fit-3.5*Mbar_err):(Mbar_err/10):(sigma_fit+3.5*Mbar_err);

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
hmle=plot(Mbar_fit,sigma_fit,'xk','MarkerSize',25);

xlabel('Mbar','FontSize',20)
ylabel('\sigma','FontSize',20)
title('Likelihood Contours and 95%, 68% Confidence Regions','FontSize',18)

xlim([Mbar_fit - 3.25*0.015,Mbar_fit + 3.25*0.015])
ylim([sigma_fit - 3.25*0.015,sigma_fit + 3.25*0.015])
