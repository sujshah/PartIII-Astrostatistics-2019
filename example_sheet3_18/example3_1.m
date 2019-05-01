clear all
close all

% optimisation options
options = optimset('MaxFunEvals',10000,'Algorithm','interior-point','Display','iter');

% read data
disp('Reading data..')
disp(' ')

fid = fopen('lensed_quasar_data.txt');

table = textscan(fid,'%f %f %f %f %f','CommentStyle','#');
fclose(fid);

ts = table{1};
y1s = table{2};
y2s = table{4};
nobs = length(ts);

y1errs = table{3};
y2errs = table{5};

figure(1)
errorbar(ts, y1s,y1errs,'.k','MarkerSize',20)
hold on
errorbar(ts, y2s,y2errs,'.b','MarkerSize',20)
hold off
title('Doubly Imaged Quasar Time Series')
xlabel('Days')
ylabel('Magnitude')
set(gca,'YDir','Reverse')
legend({'y1(t)','y2(t)'},'Location','NorthEast')

%% fit OU process to y1 data

disp('Fitting OU to y1')
disp(' ')

obj_ou1 = @(pp) -1*loglkhd_ou_process(y1s,y1errs,ts,pp(1),pp(2),pp(3));

start = [-100,10,1000]';
obj_ou1(start);

lb = [-Inf,1e-6,1e-6]';
ub = Inf(3,1);

[out,fval,exitflag,output,lambda,grad,hessian] = fmincon(obj_ou1,start,[],[],[],[],lb,ub,[],options);

minval = fval
fit_ou_pars1 = out;

c1 = fit_ou_pars1(1)
A1 = fit_ou_pars1(2)
tau1 = fit_ou_pars1(3)

tgrid = (100:1:1000)';
T_grid_mat = tgrid(:,ones(1,length(tgrid)));

% plot GP prior draws vs. dataset y1
figure(2)
for i=1:3
  rand_lc = mvnrnd(tgrid*0+c1, A1^2 *exp(-abs(T_grid_mat-T_grid_mat')/tau1));
  plot(tgrid,rand_lc,'LineWidth',2)
  hold on
end

errorbar(ts, y1s,y1errs,'.k','MarkerSize',20)
xlabel('Days')
ylabel('Magnitude')
set(gca,'YDir','Reverse')
title('Brightness Time Series y1(t) & GP Prior Draws')
hold off

%% compute posterior inference of GP curve for y1

disp('Compute GP posterior fit to y1')
disp(' ')

[condE,condCov] = gp_predict_ou(y1s,ts,y1errs,tgrid,c1,A1,tau1);

% plot dataset y1
figure(3)
errorbar(ts, y1s,y1errs,'.k','MarkerSize',20)
xlabel('Days')
ylabel('Magnitude')
set(gca,'YDir','Reverse')
title('Brightness Time Series y1(t) & GP Fit')
hold on

%plot posterior GP fit
plot(tgrid,condE,'-k','LineWidth',2)
[tvs,yvs] = errsnake(tgrid,[condE+sqrt(diag(condCov)),condE-sqrt(diag(condCov))]);
fill(tvs,yvs,[0.,0.5,0.5],'EdgeColor','none','FaceAlpha',0.5);
hold off

%% same thing for Y2

disp('Fitting OU to y2')
disp(' ')

obj_ou2 = @(pp) -1*loglkhd_ou_process(y2s,y2errs,ts,pp(1),pp(2),pp(3));
obj_ou2(start);

[out,fval,exitflag,output,lambda,grad,hessian] = fmincon(obj_ou2,start,[],[],[],[],lb,ub,[],options);

minval = fval
fit_ou_pars2 = out;

c2 = fit_ou_pars2(1)
A2 = fit_ou_pars2(2)
tau2 = fit_ou_pars2(3)

tgrid = (100:1:1000)';
T_grid_mat = tgrid(:,ones(1,length(tgrid)));

% plot GP prior draws vs. dataset y2
figure(4)
for i=1:3
  rand_lc = mvnrnd(tgrid*0+c2, A2^2 *exp(-abs(T_grid_mat-T_grid_mat')/tau2));
  plot(tgrid,rand_lc,'LineWidth',2)
  hold on
end

errorbar(ts, y2s,y2errs,'.k','MarkerSize',20)
xlabel('Days')
ylabel('Magnitude')
set(gca,'YDir','Reverse')
title('Brightness Time Series y2(t) & GP Prior Draws')
hold off

%% compute posterior inference of GP curve for y2
disp('Compute GP posterior fit to y2')
disp(' ')

[condE,condCov] = gp_predict_ou(y2s,ts,y2errs,tgrid,c2,A2,tau2);

% plot dataset y2
figure(5)
errorbar(ts, y2s, y2errs,'.k','MarkerSize',20)
xlabel('Days')
ylabel('Magnitude')
set(gca,'YDir','Reverse')
title('Brightness Time Series y2(t) & GP Fit')
hold on

%plot posterior GP fit
plot(tgrid,condE,'-k','LineWidth',2)
[tvs,yvs] = errsnake(tgrid,[condE+sqrt(diag(condCov)),condE-sqrt(diag(condCov))]);
fill(tvs,yvs,[0.,0.5,0.5],'EdgeColor','none','FaceAlpha',0.5);
hold off


%% preliminary estimates of parameters

disp('Preliminary Estimates of Parameters')
disp(' ')
% estimate mean level of OU process for y1
c_est = c1

% estimate of the overall difference in magnitudes 
% btw y2 relative to y1
dm_est = c2 - c1

% estimate the amplitude and timescale by averaging
% the (consistent) estimates from y1 and y2 separately
A_est = 0.5*( A1 + A2)
tau_est = 0.5*( tau1 + tau2)

%%
%obj = @(pp) -1*loglkhd_timedelay(y1s,y1errs,y2s,y2errs,ts,pp(1),pp(2),pp(3),pp(4),pp(5));

obj = @(pp) -1*loglkhd_timedelay(y1s,y1errs,y2s,y2errs,ts,pp(1),pp(2),pp(3),A_est,tau_est);
%obj = @(pp) -1*loglkhd_timedelay(y1s,y1errs,y2s,y2errs,ts,75,0.1,pp(1),pp(2),pp(3));
%obj = @(pp) -1*loglkhd_timedelay(y1s,y1errs,y2s,y2errs,ts,pp(1),pp(2),0,0.01,50);

%start = [50,dm_est,c_est,A_est,tau_est]';
%start = [c_est,A_est,tau_est]';
start = [70,dm_est,c_est];

obj(start)

%lb = [-Inf,-Inf,-Inf,1e-6,1e-6]';
%ub = Inf(5,1);

lb = [-Inf,-Inf,-Inf]';
ub = Inf(3,1);

%lb = [-Inf,1e-6,1e-6]';
%ub = Inf(3,1);

%lb = [-Inf, -Inf];
%ub = [Inf, Inf];

[out,fval,exitflag,output,lambda,grad,hessian] = fmincon(obj,start,[],[],[],[],lb,ub,[],options);

minval = fval;
fitpars = out;

dt_est = fitpars(1)
dm_est = fitpars(2)
c_est = fitpars(3)


%% evaluate GP marginal likelihood as fcn of dt, dm
dm_grid = (0.085:0.001:0.115)';
dt_grid = (71.5:0.01:72.5)';

[X,Y] = meshgrid(dt_grid,dm_grid);

loglkhd_grid = X*0;

for i=1:length(dm_grid)
    for j=1:length(dt_grid)
        
        loglkhd_grid(i,j) = -1*obj([dt_grid(j),dm_grid(i),c_est]);
    end
end

loglkhd_grid = loglkhd_grid - max(max(loglkhd_grid));

%% make likelihood contours
figure(6)
contour(X,Y,exp(loglkhd_grid),[0.01,0.05:0.1:0.95,0.99],'LineWidth',2)
xlabel('Time Delay \Deltat')
ylabel('Mag Shift \Deltam')
title('Marginal Likelihood Contours')

%% run Metropolis-within-Gibbs MCMC

n_mc = 1e4;

acc = zeros(3,1);

mc = zeros(n_mc,3);

disp('Initialising MCMC ...')
disp(' ')

% initialise chain
start = [73,dm_est,c_est]';
theta = start
log_post_curr = -1*obj(theta)

mc(1,:) = theta;

% proposal jump scale in each dimension (dt, dm, c)
jumps = [0.75; 0.01; 0.15];

tic;

for i=2:n_mc
   
   if mod(i,n_mc/10)==0
       disp(['i = ' num2str(i,'%.0f') ' : logpost = ' num2str(log_post_curr) ' : dt = ' num2str(theta(1))])
   end
   for j=1:length(jumps)
       
       theta_prop = theta;
       theta_prop(j) = theta(j) + jumps(j)*randn;
       
       log_post_prop = -1*obj(theta_prop);
       logr = log_post_prop - log_post_curr;
       
       if log(rand) < logr
           theta = theta_prop;
           acc(j) = acc(j)+1;
           log_post_curr = log_post_prop;
       else 
           % theta stays the same
       end
   end
   mc(i,:) = theta;
    
end

acc = acc/n_mc

runtime = toc

%% load pre-run MCMC chains (4)
%analyse chains

load mcmc_all.mat


%%
mc = mcmc_combine(mc_all,0.2*n_mc);

dt_est = mean(mc(:,1))
dm_est = mean(mc(:,2))
c_est = mean(mc(:,3))

%% plot deshifted lcs

figure(7)
errorbar(ts, y1s,y1errs,'.k','MarkerSize',20)
hold on
errorbar(ts-dt_est, y2s-dm_est,y2errs,'.b','MarkerSize',20)
hold off
title('Doubly Imaged Quasar Time Series - deshifted')
xlabel('Days')
ylabel('Magnitude')
set(gca,'YDir','Reverse')
legend({'y1(t)','y2(t-\Deltat) - \Deltam'},'Location','NorthEast')

%% posterior inference of latent light curve using combined time series

ys_comb = [y1s; y2s-dm_est];
ts_comb = [ts; ts-dt_est];
yerrs_comb = [y1errs; y2errs];

[condE,condCov] = gp_predict_ou(ys_comb,ts_comb,yerrs_comb,tgrid,c_est,A_est,tau_est);

hold on
plot(tgrid,condE,'-k','LineWidth',2)
[tvs,yvs] = errsnake(tgrid,[condE+sqrt(diag(condCov)),condE-sqrt(diag(condCov))]);
fill(tvs,yvs,[0.,0.5,0.5],'EdgeColor','none','FaceAlpha',0.5);
hold off

