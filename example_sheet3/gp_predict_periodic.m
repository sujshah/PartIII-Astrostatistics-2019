function [condE,condCov] = gp_predict_periodic(ys,ts,W,tgrid,mu,Hyps)

a = Hyps(1);
l_tilde = Hyps(2);
period = Hyps(3);

tvec = [ts;tgrid];

n_obs = length(ts);

n_tot = length(tvec);

n_grid = length(tgrid);

qobs = (1:n_obs)';

qgrid = ((n_obs+1):n_tot)';

t_mat = tvec(:,ones(1,n_tot));

CMat = a^2 *exp(-2/l_tilde^2 *sin(pi*(t_mat-t_mat')/period).^2);

Sigma = CMat;
%Sigma(qobs,qobs) = Sigma(qobs,qobs) + W + sig2*eye(n_obs);
Sigma(qobs,qobs) = Sigma(qobs,qobs) + W;

Sigma_oo = Sigma(qobs,qobs);
Sigma_og = Sigma(qobs,qgrid);
Sigma_go = Sigma(qgrid,qobs);
Sigma_gg = Sigma(qgrid,qgrid);
invSigma_oo = inv(Sigma_oo);

condE = mu*ones(n_grid,1) + Sigma_go*invSigma_oo*(ys - mu*ones(n_obs,1));

condCov = Sigma_gg - Sigma_go*invSigma_oo*Sigma_og;

end