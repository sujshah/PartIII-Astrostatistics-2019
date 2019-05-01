function Rhat = mcmc_calcrhat(mc)

s = size(mc);
nsmpl = s(1);
npars = s(2);
nchain = s(3);

if (nchain < 2)
    error('Need more than one chain!');
end

B = zeros(npars,1);
W = zeros(npars,1);
Rhat = zeros(npars,1);


chainmeans = squeeze(mean(mc,1));

%size(chainmeans)

meanchainmeans = mean(chainmeans,2);

%size(meanchainmeans)

chainvars = squeeze(var(mc,0,1));

%size(chainvars)

B = nsmpl/(nchain-1) * sum( (chainmeans - meanchainmeans(:,ones(1,nchain))).^2,2);

%B(ipar) = nsmpl./(nchain-1.) .* sum( (chainmeans-mean(chainmeans)).^2 );
%W(ipar) = mean(chainvars);

W = mean(chainvars,2);

%Rhat(ipar) = sqrt( (nsmpl-1)./nsmpl .* W(ipar) + B(ipar)./nsmpl)./sqrt(W(ipar));

%size(B)
%size(W)

Rhat = ( W*(nsmpl-1)/nsmpl + B/nsmpl )./W;

%size(Rhat)
Rhat = sqrt(Rhat);

    
end
