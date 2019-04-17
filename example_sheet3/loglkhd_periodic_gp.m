function logl = loglkhd_periodic_gp(tobs,mobs,merr,mu,a,l_tilde,period)

Tobs_mat = tobs(:,ones(1,length(tobs)));

Wmat = diag(merr.^2);

CovMat = Wmat + a^2 *exp(-2/l_tilde^2 *sin(pi*(Tobs_mat-Tobs_mat')/period).^2);

logl = -0.5*log(2*pi*det(CovMat)) - 0.5*(mu-mobs)'*(CovMat\(mu-mobs));



end