function L=logL_regress_gx1(y,yerr,X,Xerr,pp)

n = size(y,1);
k = size(X,2);

alpha = pp(1);
betas = pp(2:(k+1));
sigma = pp(k+2);

mus = pp(k+2+(1:k));
taus = pp(2*k+2 + (1:k));

if (sigma < 0) || (prod(taus) < 0)
    L = -1e6;
    return
else
    
    sigma2 = sigma^2;
    
    Vs = zeros(n,1);
    Es = zeros(n,1);
    Sigma_xs = zeros(k,k);
    
    T_plus_Sigma_x = zeros(k,k);
    
    T = diag(taus.^2);
    L = 0;
    
    for s=1:n
        
        xs = X(s,:)';
        
        Sigma_xs = diag(Xerr(s,:)'.^2);
        
        T_plus_Sigma_x = T + Sigma_xs;
        
        Es(s) = alpha + betas'*mus + (betas'*T)*( T_plus_Sigma_x \ (xs-mus) );
        
        Vs(s) = betas'*T*betas + sigma2 + yerr(s)^2 - (betas'*T)*( T_plus_Sigma_x  \ (T*betas));
        
        % add in part from likelihood of x
        L = L -0.5*( xs-mus )'*( T_plus_Sigma_x \ (xs-mus) ) - 0.5*log(det(2*pi*T_plus_Sigma_x));
        
        % add in part from likelihood of y | x
        L = L -0.5*( y(s)-Es(s) ).^2 /Vs(s) -0.5*log(2*pi*Vs(s));
        
    end
    
    
end
