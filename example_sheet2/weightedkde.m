function [kde_pdf, mean_x, mode_x, std_x, hpd68,hpd90,pactuals,hpd_lvls] = weightedkde(xs,lnw,xgrid)

%normalize weights just in case they aren't
% if they are already normalized, this doesn't change them.
lnw_norm = lnw - log(sum(exp(lnw)));

% compute effective sample size
ess = (sum(exp(lnw))^2)/sum((exp(lnw).^2));

mean_x = sum(xs .* exp(lnw_norm));

std_x = sqrt(sum( (xs-mean_x).^2 .* exp(lnw_norm)));

bw_opt = 1.06 * std_x / ess^(1/5);

n_xs = length(xs);
n_grid = length(xgrid);

density_matrix = nan(n_xs,n_grid);

for i =1:n_xs
    
    % give each host a Gaussian centered at x with std.dev given by
    % the bandwidth and then weight by the importance
    density_matrix(i,:) = exp(lnw_norm(i)) * normpdf(xgrid,xs(i),bw_opt);
    
end

%assume a regular grid, this is the grid spacing
hx = xgrid(2)-xgrid(1);

% the kde estiamte is the sum of this matrix over the hosts;
kde_pdf = sum(density_matrix,1)';
kde_mean = sum(xgrid .* kde_pdf * hx);
kde_std = sqrt(sum(hx * kde_pdf .* (xgrid-kde_mean).^2));

%check that kde integrates to one.
should_be_one = sum(kde_pdf * hx);

if (abs(should_be_one - 1) > 0.02)
    should_be_one
    error('Unitarity problem')
end

% find perc% credible interval
%hpd68 = 0;
%perc = 0.68;

maxpdf = max(kde_pdf);

levels = (0.1:0.01:0.9)'*maxpdf;
ps = levels*0;

qmax = find(kde_pdf==maxpdf);

left_pdf = kde_pdf(1:qmax);
left_x = xgrid(1:qmax);

% find values of x corresponding to density levels on the left side
lefts = interp1(left_pdf,left_x,levels);

right_pdf = kde_pdf(qmax:end);
right_x = xgrid(qmax:end);

% find values of x corresponding to density levels on the right side
rights = interp1(right_pdf,right_x,levels);

for i=1:length(levels);
    % find the samples in this range defined by the level set
   xq = (xs > lefts(i) & xs < rights(i));
   % calculate
   ps(i) = sum( exp(lnw_norm(xq)));
end

%[lefts, rights, ps]

% figure
% plot(ps,levels/maxpdf)
% xlabel('posterior probability enclosed in HPD interval')
% ylabel('density level / max pdf')

percs = [0.68, 0.90];

n_percs = length(percs);
hpds = nan(n_percs,2);
pactuals = nan(n_percs,1);
hpd_lvls = nan(n_percs,1);

for j=1:length(percs)
    % desired hpd percent interval
    perc = percs(j);
    best = find( abs(ps-perc) == min(abs(ps-perc)) );
    best = best(1);
    pactuals(j) = ps(best);
    hpds(j,:) = [lefts(best), rights(best)];
    hpd_lvls(j) = levels(best);
end

mode_x = xgrid(kde_pdf == max(kde_pdf));

hpd68 = hpds(1,:);
hpd90 = hpds(2,:);

% makes a plot
% to turn off plot, set doplot=0
doplot = 1
if doplot == 1
    
    plot(xgrid,kde_pdf,'-k','LineWidth',2)
    
    hold on
    
    % plot mean
    hmn=plot([ mean_x,  mean_x],[0,1.1*max(kde_pdf)],'--k','LineWidth',2);
    
    % plot mode
    hmd = plot([ mode_x,  mode_x],[0,1.1*max(kde_pdf)],'-k','LineWidth',2);
    
    % plot 68% hpd interval
    h68=plot([hpd68(1), hpd68(1)],[0, hpd_lvls(1)],'--r','LineWidth',2);
    plot([hpd68(2), hpd68(2)],[0, hpd_lvls(1)],'--r','LineWidth',2)
    
    % plot 90% hpd interval
    h90=plot([hpd90(1), hpd90(1)],[0, hpd_lvls(2)],'-r','LineWidth',2);
    plot([hpd90(2), hpd90(2)],[0, hpd_lvls(2)],'-r','LineWidth',2)
    
    
    legend([hmn,hmd,h68,h90],...
        {['Mean = ' num2str(mean_x,'%.2f')],...
        ['Mode = ' num2str(mode_x,'%.2f')],...
        ['68% hpd = [ ' num2str(hpd68(1),'%.2f') ', ' num2str(hpd68(2),'%.2f') ' ]'],...
        ['90% hpd = [ ' num2str(hpd90(1),'%.2f') ', ' num2str(hpd90(2),'%.2f') ' ]']},...
        'Location','NorthEast','FontSize',16,'Box','Off')
    
    xlabel('X','FontSize',18)
    ylabel('posterior pdf','FontSize',18)
    hold off
end


end