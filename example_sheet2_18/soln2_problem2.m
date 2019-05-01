clear all
close all

% set observed data of the Large Magellanic Cloud (LMC)
% Patel et al. 2017, Table 1
vmax_obs = 85;
vmax_err = 10;

r_obs = 50;
r_err = 5;

vtot_obs = 321;
vtot_err = 24;

j_obs = 15688;
j_err = 1788;

% read in prior samples of host-satellite systems from Illustris simulation
fid = fopen('Patel17b_Illustris_Data_KM.txt');

% columns 2,4,6,7,8 are:
% host virial mass, massive satellite vmax, position, velocity, and magnitude of j

table = textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f','CommentStyle','#');
fclose(fid);

% Mass in 10^10 M_solar
Mass = table{2};

N = length(Mass)

%calculate log10(Mass/M_solar)
log10Mass = 10 + log10(Mass);

% vmax = rotational velocity of satellite
vmaxs = table{4};

% position = distance from central galaxy
rs = table{6};

% vtot = velocity of satellite
vtots = table{7};

% j = angular momentum of satellite
js = table{8};

%%

% calculate log likelihoods for vmax, r, vtot for all prior
% samples

loglkhd_vmaxs = -0.5*(vmax_obs - vmaxs).^2 ./vmax_err^2 -0.5*log(2*pi*vmax_err^2);
loglkhd_rs    =  -0.5*(r_obs - rs).^2 ./ r_err^2 - 0.5*log(2*pi*r_err^2) ;
loglkhd_vtots = -0.5*(vtot_obs - vtots).^2 ./ vtot_err^2 - 0.5*log(2*pi*vtot_err^2);
loglkhd_js = -0.5*(j_obs - js).^2 ./ j_err^2 -0.5*log(2*pi*vtot_err^2);

% calculate combined log likelihood for (vmax, r, vtot)
loglkhds = loglkhd_vmaxs + loglkhd_js
%loglkhds = loglkhd_vmaxs + loglkhd_rs + loglkhd_vtots;

%ln weights unnormalized.  Importance weights are propto likelihoods
lnw = loglkhds;

%ln weights normalized: ensures sum of weights equals 1
lnw_norm = lnw - log(sum(exp(lnw)));

% calculate effective sample size
ess = (sum(exp(lnw_norm))^2)/sum((exp(lnw_norm).^2))

figure(1)
histogram(lnw_norm(exp(lnw_norm) > 1e-1000))
xlabel('log Importance Weights','FontSize',16)
ylabel('Number','FontSize',16);
title('Distribution of Importance Weights','FontSize',16)
text(-500,1200,['Prior samples N = ' num2str(N,'%.0f')],'FontSize',16)
text(-500,1000,['Effective SS = ' num2str(ess,'%.2f')],'FontSize',16)

%% calculate posterior in log10 Mass
logMass_grid = (10:0.01:14.5)' ;

% plot prior KDE
figure(2)
[prior_kde, log10Mi] = ksdensity(log10Mass);

plot(log10Mi,prior_kde,'-b','LineWidth',2)
hold on

[kde_pdf_lM, mean_lM, mode_lM, std_lM, hpd68_lM, hpd90_lM, p_actual_lM, hpd_lvls_lM] = weightedkde(log10Mass,lnw,logMass_grid);

mean_lM
mode_lM
std_lM
hpd68_lM
hpd90_lM
p_actual_lM

xlim([10.5,14.5])
hold off
xlabel('log_{10}[ M_{vir} / M_{solar} ]','FontSize',12)
ylabel('Probability Density','FontSize',12)

text(11.2,0.75,'Prior','FontSize',12,'Color','b')
text(11.5,2,'Posterior','FontSize',12) 
title('Bayesian Inference of MW Mass','FontSize',12)

%% calculate in linear Mass scale
% Mass_grid = (0:1:1500)';
% 
% figure(3)
% [kde_pdf_M, mean_M, mode_M, std_M, hpd68_M, hpd90_M, p_actual_M, hpd_lvls_M] = weightedkde(Mass,lnw,Mass_grid);
% xlabel('M_{vir} (10^{10} M_{solar})','FontSize',18)
% mean_M
% mode_M
% std_M
% hpd68_M
% w68_M = hpd68_M(2) - hpd68_M(1)
% hpd90_M
% p_actual_M
