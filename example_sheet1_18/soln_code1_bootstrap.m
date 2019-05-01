clear all
close all;

% number of bootstraps to run
nboot = 1000;

filename = 'example_sheet1_prb1_data.txt';

disp(' ')
disp(['Reading in data from ' filename])

fid = fopen(filename);

table = textscan(fid,'%f %f %f %f %f','CommentStyle','#');

% unpack the data
N = length(table{1});
mhats = table{2};
sigma_ms = table{3};

muhats = table{4};
sigma_mus = table{5};

Mobs_orig = mhats - muhats;
sigmas_orig = sqrt(sigma_mus.^2 + sigma_ms.^2);

% create bootstrap samples
Mobs_boot = zeros(N,nboot);
sigmas_boot = zeros(N,nboot);

[~, bootsam] = bootstrp(nboot,@(x) x, Mobs_orig);

for i=1:nboot
    Mobs_boot(:,i) = Mobs_orig(bootsam(:,i));
    sigmas_boot(:,i) = sigmas_orig(bootsam(:,i));
end
    
% plot the data
if nboot < 5
    
    figure(1);
    
    for i=1:nboot
        
        subplot(2,2,i)
        
        edges = -18.7:0.1:-17;
        
        hobs=histogram(Mobs_boot(:,i),edges);
        set(gca,'FontSize',18)
        hold on;
        
        xlabel('Absolute Magnitude','FontSize',20)
        ylabel('Number','FontSize',20)
        %legend([hobs],{'Observed Abs Mag'},'FontSize',18,'Location','NorthEast')
        title(['N = ' num2str(N,'%.0f') ', Sample STD(obs) = ' num2str(std(Mobs_boot(:,i)),'%.3f')],'FontSize',18)
        
        xlim([min(edges),max(edges)])
        ymax = 1+max(hobs.Values);
        ylim([0,ymax])
        
    end
    
end

%% fit data

mu_boots = zeros(nboot,1);
sigma_int_boots = zeros(nboot,1);

disp(' ')
disp('Finding Maximum Likelihood Solution')

for i=1:nboot
    Mobs = Mobs_boot(:,i);
    sigmas = sigmas_boot(:,i);
    
    % a quick function of theta = (mu, sigma_int) to evaluate the negative log likelihood
    neglogL = @(theta) 0.5* sum( (Mobs-theta(1)).^2 ./(sigmas.^2 + theta(2)^2) + log(2*pi*(sigmas.^2 + theta(2)^2)) );
    
    % set some options
    options=optimset('Display','None','UseParallel',0,'Algorithm','active-set');
    
    init = [mean(Mobs),std(Mobs)];
    
    lb = [-Inf;0];
    ub = [Inf;Inf];
    
    % optimise the function to get the MLE
    % also outputs the hessian of the negative log lkhd at the MLE
    % this matrix is the observed fisher information
    [bestfit,fval,exitflag,output,lambda,grad,hessian]= fmincon(neglogL,init,[],[],[],[],lb,ub,[],options);
    
    Mbar_fit = bestfit(1);
    sigma_fit = bestfit(2);
    
    % invert the fisher matrix to get an asymptotic approximation to the
    % covariance matrix
    cov_unc = inv(hessian);
    
    % the diagonals are lower bounds on the variance of the MLEs
    param_errs = sqrt(diag(cov_unc));
    
    Mbar_err = param_errs(1);
    sigma_err = param_errs(2);
    
    mu_boots(i) = Mbar_fit;
    sigma_int_boots(i) = sigma_fit;
    
end

figure(2)
subplot(2,1,1)
histogram(mu_boots)
set(gca,'FontSize',18)
title([num2str(nboot,'%.0f') ' bootstraps, Mbar_{hat} = ' num2str(mean(mu_boots),'%.3f') ' \pm ' num2str(std(mu_boots),'%.3f')]) 

subplot(2,1,2)
histogram(sigma_int_boots)
xx = sigma_int_boots > 0.00;
set(gca,'FontSize',18)
title(['\sigma_{int,hat} = ' num2str(mean(sigma_int_boots(xx)),'%.3f') ' \pm ' num2str(std(sigma_int_boots(xx)),'%.3f')]) 

% disp(' ')
% disp('-----Maximum Likelihood Results-----')
% disp(['mu_hat = ' num2str(mu_fit,'%.3f') ' +/- ' num2str(mu_err,'%.3f') ' mag'])
% disp(['sigma_hat = ' num2str(sigma_fit,'%.3f') ' +/- ' num2str(sigma_err,'%.3f') ' mag'])
% disp('------------------------------------')
% disp(' ')

