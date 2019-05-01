function l = loglkhd_ou_process(y1s,y1errs,tobs,c,A,tau)

    Ys = y1s;
    ts = tobs;
    one = Ys*0 + 1;
    
    W = diag(y1errs.^2);
    
    T_grid_mat = ts(:,ones(1,length(ts)));
    
    totalcov = A^2 *exp(-abs(T_grid_mat-T_grid_mat') /tau) + W;
    
    %detcov = det(totalcov)
    %invcov = inv(totalcov);
    cholL = chol(totalcov,'lower');
    logdetcov = 2 * sum(log(diag(cholL)));
    
    resid = (Ys - one*c);
    
    l = -0.5*logdetcov - 0.5*resid' * (cholL' \(cholL\(resid)));
    %l = -0.5*logdetcov - 0.5*resid'*inv(totalcov)*resid;
    


end