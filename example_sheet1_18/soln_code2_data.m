clear all
close all;

filename = 'example_sheet1_prb2_data.txt';

disp(['Reading in data from ' filename])

fid = fopen(filename);

table = textscan(fid,'%s %s %f %f','CommentStyle','#');

% unpack the data
N = length(table{1});
sns = table{1};
Ohats = table{3};
sigmas = table{4};

% plot the data
figure(1);

hobs=histogram(Ohats,10);
set(gca,'FontSize',20)

xlabel('Observed Apparent Colour','FontSize',20)
ylabel('Number','FontSize',20)

title(['N = ' num2str(N,'%.0f') ' SNe Ia from Jha, Riess & Kirshner 2007'],'FontSize',18)


%% fit data

disp('Finding Maximum Likelihood Solution')

% a quick function of theta = (mu, sigma_int, tau) to evaluate the negative log likelihood 
neglogL = @(theta) -1*sum( -log(theta(3)) ...
    + 0.5*( sqrt(sigmas.^2+theta(2)^2)/theta(3)).^2 ...
    - (Ohats-theta(1))/theta(3) ...
    + log(normcdf( (Ohats-theta(1))./sqrt(sigmas.^2+theta(2)^2) - sqrt(sigmas.^2+theta(2)^2)/theta(3))));


% set some options
options=optimset('Display','iter','UseParallel',0,'Algorithm','active-set');
   
init = [mean(Ohats),std(Ohats),0.4];

lb = [-Inf;0;0.01];
ub = [Inf;Inf;Inf];

% optimise the function to get the MLE
% also outputs the hessian of the negative log lkhd at the MLE 
% this matrix is the observed fisher information
[bestfit,fval,exitflag,output,lambda,grad,hessian]= fmincon(neglogL,init,[],[],[],[],lb,ub,[],options);

mu_fit = bestfit(1)
sigma_fit = bestfit(2)
tau_fit = bestfit(3)

% invert the fisher matrix to get an asymptotic approximation to the
% covariance matrix
cov_unc = inv(hessian)

% the diagonals are lower bounds on the variance of the MLEs
param_errs = sqrt(diag(cov_unc));

mu_err = param_errs(1)
sigma_err = param_errs(2)
tau_err = param_errs(3)

disp(' ')
disp('-----Maximum Likelihood Results-----')
disp(['mu_hat = ' num2str(mu_fit,'%.3f') ' +/- ' num2str(mu_err,'%.3f') ' mag'])
disp(['sigma_int_hat = ' num2str(sigma_fit,'%.3f') ' +/- ' num2str(sigma_err,'%.3f') ' mag'])
disp(['tau_hat = ' num2str(tau_fit,'%.3f') ' +/- ' num2str(tau_err,'%.3f') ' mag'])
disp('------------------------------------')
disp(' ')

%%
% plot bestfit solution

figure(2)
O_grid = (0.8:0.01:2.0)';

% evaluates sampling density using MLE parameters
neglogL1 = @(O) -1*( -log(bestfit(3)) ...
    + 0.5*( sqrt(bestfit(2)^2)/bestfit(3)).^2 ...
    - (O-bestfit(1))/bestfit(3) ...
    + log(normcdf( (O-bestfit(1))./sqrt(bestfit(2)^2) - sqrt(bestfit(2)^2)/bestfit(3))));

hobs=histogram(Ohats,10,'Normalization','pdf');
set(gca,'FontSize',20)
box on
hold on

xlabel('Observed Apparent Colour','FontSize',20)
ylabel('Number','FontSize',20)

title(['N = ' num2str(N,'%.0f') ' SNe Ia'],'FontSize',18)
hmod=plot(O_grid,exp(-1*neglogL1(O_grid)),'-r','LineWidth',2);

legend([hobs,hmod],{'Observed Colours','P(O_{hat} | \mu_C, \sigma_C, \tau)'},'Box','Off')

text(1.3,3.5,['\mu_{hat} = ' num2str(mu_fit,'%.3f') ' \pm ' num2str(mu_err,'%.3f') ' mag'],'FontSize',20)
text(1.3,3.0,['\sigma_{int,hat} = ' num2str(sigma_fit,'%.3f') ' \pm ' num2str(sigma_err,'%.3f') ' mag'],'FontSize',20)
text(1.3,2.5,['\tau_{hat} = ' num2str(tau_fit,'%.3f') ' \pm ' num2str(tau_err,'%.3f') ' mag'],'FontSize',20)
hold off

