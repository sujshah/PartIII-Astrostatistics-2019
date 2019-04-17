close all
clear all

fs = 12;

% read data
fid = fopen('variable_star.txt');

tab = textscan(fid,'%f %f %f','CommentStyle','#');

tobs = tab{1};
mobs = tab{2};
merr = tab{3};

figure(1)
errorbar(tobs,mobs,merr,'.','MarkerSize',fs)
xlabel('Time (days)','FontSize',fs);
ylabel('Magnitude','FontSize',fs);
set(gca,'FontSize',fs);
set(gca,'YDir','Reverse')

%% assume a zero mean GP
mu = 0;

%% compute profile likelihood as a function of T

trial_Ts = 10:1:1000;
%trial_Ts = 190:0.1:210;
n_trials = length(trial_Ts);
logLs = zeros(n_trials,1);
params = zeros(n_trials,2);
perrs = zeros(n_trials,2);

for i=1:length(trial_Ts)
    if mod(i,n_trials/10)==0
       disp(['Trial period i = ' num2str(i)]) 
    end
   
    obj = @(pp) -1*loglkhd_periodic_gp(tobs,mobs,merr,mu,pp(1),pp(2),trial_Ts(i));

    start = [1,1];
    low_bnds = [0, 0];
    upp_bnds = [100,100];

    [out,fval,exitflag,output,lambda,grad,hessian] = fmincon(obj,start,[],[],[],[],low_bnds,upp_bnds,[]);

    logLs(i) = -fval;
    params(i,:) = out;
    fisher = inv(hessian);
    perrs(i,:) = sqrt(diag(fisher));

end

%%

figure(2)
plot(trial_Ts,logLs);
xlabel('Trial Period T','FontSize',fs);
ylabel('Log Profile Likelihood','FontSize',fs);

figure(3)
plot(trial_Ts,exp(logLs - max(logLs)),'LineWidth',2);
xlabel('Trial Period T','FontSize',fs);
ylabel('Profile Likelihood','FontSize',fs);

%% optimize full log likelihood

obj = @(pp) -1*loglkhd_periodic_gp(tobs,mobs,merr,mu,pp(1),pp(2),pp(3));

besti = find(logLs == max(logLs));

% start optimisation in 3 parameters from the best values at the 
% peak of the profile likelihood
start = [params(besti,:), trial_Ts(besti)];
low_bnds = [0, 0, 0];
upp_bnds = [100,100,1000];

[out,fval,exitflag,output,lambda,grad,hessian] = fmincon(obj,start,[],[],[],[],low_bnds,upp_bnds,[]);

%% compute maximum likelihood and Fisher matrix

logL_best = -fval;
param_best = out
fisher = inv(hessian);
perrs_best = sqrt(diag(fisher))'

%% compute interpolated and future values;

tgrid = (1:2000)';

[condE,condCov] = gp_predict_periodic(mobs,tobs,diag(merr.^2),tgrid,mu,out);

figure(4)
errorbar(tobs,mobs,merr,'.','MarkerSize',fs)
xlabel('Time (days)','FontSize',fs);
ylabel('Magnitude','FontSize',fs);
set(gca,'FontSize',fs);
set(gca,'YDir','Reverse')
hold on
plot(tgrid,condE)
condStd = sqrt(diag(condCov));
hold on
[tvs,yvs] = errsnake(tgrid,[condE+condStd,condE-condStd]);
fill(tvs,yvs,[0.,0.5,0.5],'EdgeColor','none','FaceAlpha',0.5);
hold off

condE(1800)
condStd(1800)

%% plot phase folde data vs. GP

figure(5)
errorbar(mod(tobs,param_best(3)),mobs,merr,'.','MarkerSize',fs)
xlabel('Phase = Folded Time [t_{obs} mod T_{mle} ] (days)','FontSize',fs);
ylabel('Magnitude','FontSize',fs);
set(gca,'FontSize',fs);
set(gca,'YDir','Reverse')
hold on
plot(tgrid,condE)
hold on
[tvs,yvs] = errsnake(tgrid,[condE+condStd,condE-condStd]);
fill(tvs,yvs,[0.,0.5,0.5],'EdgeColor','none','FaceAlpha',0.5);
xlim([0,200])
hold off

%% plot vs true lc

lc_file = 'periodic_gp_2_ok4.txt';

fid = fopen(lc_file);
tab = textscan(fid,'%f %f','CommentStyle','#');

ts = tab{1};
ms = tab{2};

figure(6)
hobs=errorbar(mod(tobs,param_best(3)),mobs,merr,'.','MarkerSize',fs);
xlabel('Phase = Folded Time [t_{obs} mod T_{mle} ] (days)','FontSize',fs);
ylabel('Magnitude','FontSize',fs);
set(gca,'FontSize',fs);
set(gca,'YDir','Reverse')
hold on
hgp=plot(tgrid,condE);
hold on
[tvs,yvs] = errsnake(tgrid,[condE+condStd,condE-condStd]);
fill(tvs,yvs,[0.,0.5,0.5],'EdgeColor','none','FaceAlpha',0.5);
xlim([0,200])

htrue=plot(ts,ms-1.1805,'-k','LineWidth',3);
hold off

legend([hobs,hgp,htrue],{'Obs Data','Post GP','True'},'Location','NorthWest')

