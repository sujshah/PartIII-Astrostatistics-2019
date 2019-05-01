function [condE,condCov] = gp_predict_ou(ys,ts,yerrs,tgrid,c,A,tau)


% if length(Hyps) == 2
% 
%     sig2 = 1;
% elseif length(Hyps) == 3
% 
%     sig2 = Hyps(3);
% end

tvec = [ts;tgrid];

n_obs = length(ts);

n_tot = length(tvec);

n_grid = length(tgrid);

qobs = (1:n_obs)';

qgrid = ((n_obs+1):n_tot)';

t_mat = tvec(:,ones(1,n_tot));

CMat = A^2 *exp(-abs(t_mat-t_mat') /tau );

W = diag(yerrs.^2);

Sigma = CMat;
%Sigma(qobs,qobs) = Sigma(qobs,qobs) + W + sig2*eye(n_obs);
Sigma(qobs,qobs) = Sigma(qobs,qobs) + W;

Sigma_oo = Sigma(qobs,qobs);
Sigma_og = Sigma(qobs,qgrid);
Sigma_go = Sigma(qgrid,qobs);
Sigma_gg = Sigma(qgrid,qgrid);
invSigma_oo = inv(Sigma_oo);

condE = c*ones(n_grid,1) + Sigma_go*invSigma_oo*(ys - c*ones(n_obs,1));

condCov = Sigma_gg - Sigma_go*invSigma_oo*Sigma_og;

end