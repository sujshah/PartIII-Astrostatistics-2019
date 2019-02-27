% metropolis within gibbs sampling

clear all
close all

% fontsize
fs = 12

% set data
y1 = 0;
y2 = 0;

ydata = [y1;y2];

% known covariance matrix, unknown means theta = (theta1, theta2)
rho = 0.5;
Sigma = [1, rho; rho, 1];

% define the posterior density
posterior = @(x) mvnpdf([x(1),x(2)],ydata',Sigma);

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

%% run Metropolis-within-Gibbs sampler

% number of iterations
n_mc = 10000;

mc = zeros(n_mc,2);

% set scale of metropolis proposal dist'n
taus = [1;1];

% record acceptances of moves in each direction
acceptances = [0;0];

%initial position
theta = [-1.5,1.5];

mc(1,:) = theta;

for i=2:n_mc
    
    % alternate between Metropolis rules in each
    % parameter P(theta1 | theta2, y) or P(theta2 | theta1, y)
    
   if mod(i,2) == 0
         % Metropolis rule on theta1
         theta1_prop = theta(1) + taus(1)*randn;
         theta_prop = [theta1_prop,theta(2)];
         
         % evalutate ratio of joint posteriors at proposed and current
         r = posterior(theta_prop) / posterior(theta);
         % but this is really just the ratios of the conditionals
         % r = P(theta1_prop | theta_2, y) / P(theta1 | theta_2, y)
         
         if rand < r
             theta(1) = theta1_prop;
             acceptances(1) = acceptances(1) + 1;
         else
             %keep theta1
         end
         
   elseif mod(i,2) == 1
         % Metropolis rule on theta2
         theta2_prop = theta(2) + taus(2)*randn;
         theta_prop = [theta(1),theta2_prop];
         
         r = posterior(theta_prop) / posterior(theta);
         
         if rand < r
             theta(2) = theta2_prop;
             acceptances(2) = acceptances(2) + 1;
         else
             %keep theta2
         end
   end
   
   % record current position
   mc(i,:) = theta;
   
end

accept = acceptances/(n_mc/2)

plot(mc(:,1), mc(:,2),'-k','Linewidth',2)
title(['MCMC iterations = ' num2str(n_mc)])

hold off
%%
figure(2)
subplot(2,1,1)
plot(mc(:,1),'LineWidth',2)
ylabel('\theta_1')
title(['Accept Ratio ( \theta_1 ) = ' num2str(accept(1),'%.2f')])

subplot(2,1,2)
plot(mc(:,2),'LineWidth',2)
ylabel('\theta_2')
xlabel('MCMC Sample')
title(['Accept Ratio ( \theta_2 ) = ' num2str(accept(2),'%.2f')])

% cut out burn-in
mc = mc(1:n_mc*0.2,:);

%% plot joint
figure(3)
contour(X,Y,pdf_mesh,'LineWidth',2)
xlabel('\theta_1','FontSize',fs)
ylabel('\theta_2','FontSize',fs)
hold on
%plot(mc(:,1), mc(:,2),'-k')
plot(mc(:,1), mc(:,2),'.', 'Color',[0.5,0.5,0.5])
title([num2str(n_mc,'%.0f') ' MCMC Samples, target: P( \theta_1, \theta_2 | y)'])

%% plot marginal and histograms
tt = (-4:0.05:4)';

figure(4)
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

