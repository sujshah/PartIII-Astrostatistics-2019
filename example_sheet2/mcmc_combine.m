function y = mcmc_combine(mc,cut)

%cuts each chain by cut and combines them into one

s = size(mc);

nsmpl = s(1);
npars  = s(2);
nchains = s(3);

y = zeros(nchains*(nsmpl-cut),npars);

for i=1:nchains
    y( (i-1)*(nsmpl-cut)+1:i*(nsmpl-cut),:) = mc(cut+1:end,:,i);
end

