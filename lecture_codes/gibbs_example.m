% gibbs sampling

clear all
close all

% fontsize
fs = 12;

% set data
y1 = 0;
y2 = 0;

ydata = [y1;y2];

% known covariance matrix, unknown means theta = (theta1, theta2)
rho = 0.999;
Sigma = [1, rho; rho, 1];

% make a contour plot of the posterior density in theta = (theta1, theta2);
theta_grid = (-2:0.05:2)';

[X,Y] = meshgrid(theta_grid,theta_grid);
pdf_mesh = X*0;

for i=1:length(theta_grid)
    
    for j=1:length(theta_grid)
    
        pdf_mesh(i,j) = mvnpdf([theta_grid(i), theta_grid(j)],ydata',Sigma);
        
    end
    
end

figure(1)
contour(X,Y,pdf_mesh)
xlabel('\theta_1','FontSize',fs)
ylabel('\theta_2','FontSize',fs)
hold on
title('P(\theta_1, \theta_2 | y)')

%% run Gibbs sampler

% number of iterations
n_mc = 100;

mc = zeros(n_mc,2);

theta = [-1.5,1.5];

mc(1,:) = theta;

for i=2:n_mc
    
    % alternate between drawing from each conditional dist'n in each
    % parameter P(theta1 | theta2, y) or P(theta2 | theta1, y)
   if mod(i,2) == 0
         % draw new theta1 from P(theta1 | theta2, y)
        theta(1) = y1 + rho*(theta(2)- y2) + sqrt(1-rho^2) * randn;
   
   elseif mod(i,2) == 1
         % draw new theta2 from P(theta2 | theta1, y)
         theta(2) = y1 + rho*(theta(1)- y1) + sqrt(1-rho^2) * randn;
   end
   
   % record current position
   mc(i,:) = theta;
   
end

plot(mc(:,1), mc(:,2),'-k','Linewidth',2)
title(['MCMC iterations = ' num2str(n_mc)])
hold off

% cut out burn-in
mc = mc(1:n_mc*0.2,:);

%% plot joint
figure(2)
contour(X,Y,pdf_mesh,'LineWidth',2)
xlabel('\theta_1','FontSize',fs)
ylabel('\theta_2','FontSize',fs)
hold on
plot(mc(:,1), mc(:,2),'.k', 'Color',[0.5,0.5,0.5])
title([num2str(n_mc,'%.0f') ' MCMC Samples, target: P(\theta_1, \theta_2 | y)'])


%% plot marginal and histograms
tt = (-4:0.05:4)';

figure(3)
subplot(2,1,1)
histogram(mc(:,1),'Normalization','pdf')
hold on
plot(tt,normpdf(tt,y1,1),'LineWidth',2)
hold off
ylabel('P(\theta_1 | y)')
xlabel('\theta_1')
title([num2str(n_mc,'%.0f') ' MCMC Samples'])

subplot(2,1,2)
histogram(mc(:,2),'Normalization','pdf')
hold on
plot(tt,normpdf(tt,y2,1),'LineWidth',2)
hold off
ylabel('P(\theta_2 | y)')
xlabel('\theta_2')

