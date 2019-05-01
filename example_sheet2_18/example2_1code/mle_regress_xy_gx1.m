function [params,ee,L,fisher] = mle_regress_xy_gx1(y,yerr,X,Xerr,varargin)

% MLE estimate of linear regression with measurement errors in y and x
% Uses the "Gaussian structural model" for linear regression of 
% Kelly, B. 2007, ApJ, 665,1489 (Section 4.3) 
% using a k=1 uncorrelated gaussian prior on the x-covariates
% maximizes the log likelihood (Eq. 32-33)

% y is the vector of responses
% yerr is the vector of measurement errors on y
% X is the design matrix s.t.
% y = alpha + X beta + yerr
% with X having n rows and k columns
% n = number of objects 
% k = number of covariates excluding 1
%
% Xerr is a matrix the same size as X containing measurement errors on each
% X value.
% 
% params returns [alpha, beta, sigma, mu, tau]
% where alpha is the intercept
% beta is a k-vector of regression coefficients for each column of X
% sigma^2 is the intrinsic variance about the linear relation
% mu is a k-vector of the means in each x-covariate
% tau is a k-vector of the std deviations in each x-covariate
%
% ee is a vector of standard errors on each fit parameter derived from the Fisher
% matrix at the MLE
%
% L is the value of the log likelihood at the fit parameters
%
% total number of parameters is p = 3*k + 2;
%
% varargin is an optional p  vector with fixed values of the
% parameters or nan to float

n = size(y,1);
k = size(X,2);
p = 3*k+2;

% start with the GLS estimate
W = diag(yerr.^2);
b_init = regressW(y,[ones(n,1),X],W);
sigma_init = std(y-[ones(n,1),X]*b_init);

low_bnds = -Inf(p,1);
low_bnds(k+2) = 0;
low_bnds(2*k+2+(1:k)) = 0;

upp_bnds = Inf(p,1);

start = [b_init; sigma_init; mean(X)';std(X)'];

obj = @(pp) -1*logL(y,yerr,X,Xerr,pp);

if nargin < 4
    error('mle_regress_xy_gx1(y,yerr,X,Xerr,[fixed])')
elseif nargin == 4
    Aeq = [];
    beq = [];
    notfixed = true(p,1);
else
    fixed = varargin{1};
    if length(fixed) ~= p
        error(['fixed must have ' num2str(p) ' elements']);
    end
    
    q = find(~isnan(fixed));
    
    notfixed = isnan(fixed);
    
    if isempty(q)
        Aeq = [];
        beq = [];
    else
        Aeq = zeros(length(q),p);
        beq = zeros(length(q),1);
        for i=1:length(q)
            Aeq(i,q(i)) = 1;
            beq(i) = fixed(q(i));
        end

    end
   
    
end

%out = fminunc(obj,[0.2; 1.; 0.02.^2]);
options = optimset('MaxFunEvals',10000,'Algorithm','interior-point','Display','off');
%options = optimset('MaxFunEvals',10000,'Algorithm','interior-point','Display','iter');

[out,fval,exitflag,output,lambda,grad,hessian] = fmincon(obj,...
    start,...
    [],[],Aeq,beq,low_bnds,upp_bnds,[],options);

hessian = hessian(notfixed, notfixed);

params = out ;
fisher = inv(hessian);
errs = sqrt(diag(fisher));

ee = zeros(p,1);
ee(notfixed) = errs;


L = -fval;

end

function L=logL(y,yerr,X,Xerr,pp)

n = size(y,1);
k = size(X,2);

alpha = pp(1);
betas = pp(2:(k+1));
sigma = pp(k+2);

mus = pp(k+2+(1:k));
taus = pp(2*k+2 + (1:k));

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
