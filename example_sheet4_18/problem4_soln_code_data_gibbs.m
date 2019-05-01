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

% definite the observed absolute magnitude
Ds = mhats - muhats;

% Ms is the true value underlying Ds

sigmas = sqrt(sigma_mus.^2 + sigma_ms.^2);

% plot the data
figure(1);

edges = -18.7:0.1:-17;

hobs=histogram(Ds,edges);
hold on;

hold off
set(gca,'FontSize',20)

xlabel('Absolute Magnitude','FontSize',20)
ylabel('Number','FontSize',20)

legend([hobs],{'Observed Abs Mag'},'FontSize',18,'Location','NorthEast')
title(['N = ' num2str(N,'%.0f') ', Sample STD(obs) = ' num2str(std(Ds),'%.3f')],'FontSize',18)

xlim([min(edges),max(edges)])
ymax = 1+max(hobs.Values);
ylim([0,ymax])

%% fit data

disp('Running Gibbs Sampler...');

n_mc = 1000
n_chains = 4

mc = zeros(n_mc, N+2, n_chains);

for c = 1:n_chains
    
    % initialise / randomise parameters
    Ms = Ds + randn(N,1) .* sigmas;
    
    mu = mean(Ms);
    tau2 = var(Ms);
    
    for i=1:n_mc
        
        a = (N-3)/2;
        
        b = (N-1)/2  * var(Ms);
        
        tau2 = 1/ gamrnd(a, 1/b);
        
        mu = mean(Ms) + randn*sqrt(tau2/N);
        
        for s=1:N
            
            weighted_mean = (Ds(s)/sigmas(s)^2 + mu/tau2) / (1/sigmas(s)^2 + 1/tau2);
            comb_variance = 1/(1/sigmas(s)^2 + 1/tau2);
            
            Ms(s) = weighted_mean + randn*sqrt(comb_variance);
            
            mc(i,s,c) = Ms(s);
        end
        
        mc(i,end-1,c) = mu;
        mc(i,end,c) = sqrt(tau2);
        
    end
    
end

disp('Done!');

gr = mcmc_calcrhat(mc);

max_gr = max(gr)

%%  Here you should make trace plots of the multiple chains to check convergence, etc.

% check convergence, thinning, etc
% plot sample autocorrelation function of the tau parameter (the slowest to
% convergence to pick an thinning factor 
autocorr(mc(:,N+2,1),100)

%% thin by the thinning factor = 10 and remove the first 20% of the chain
mc = mc(0.2*n_mc:10:end,:,:);

%% combine chains 
mc_final = mc(:,:,1);
for i=2:n_chains
   mc_final = [mc_final; mc(:,:,2)];
end
mc = mc_final;

%%

post_mu_mean = mean(mc(:,N+1))
post_mu_std  = std(mc(:,N+1))
figure(1)
histogram(mc(:,N+1))
xlabel('\mu','FontSize',14)
ylabel('Posterior Probability Density','FontSize',14);

%%
post_tau_mean = mean(mc(:,N+2))
post_tau_std  = std(mc(:,N+2))
figure(2)
histogram(mc(:,N+2))
xlabel('\tau','FontSize',14)
ylabel('Posterior Probability Density','FontSize',14);




