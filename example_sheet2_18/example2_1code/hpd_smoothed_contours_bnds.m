function [contoursout, mode, X, Y, density] = hpd_smoothed_contours_bnds(data,xbnd,ybnd,nsmooth,varargin)

%finds the 95% and 68% highest posterior density contours of a bivarite
%data set

%uses kde2d.m kernel density estimation

if nargin > 2
    res = varargin{1};
else
    res = 64;
end

mn = min(data);
mx = max(data);

[bandwidth,density,X,Y]=kde2d(data,res,mn,mx);

peak = max(max(density));

mode = [ X(find(density==peak)), Y(find(density==peak)) ];

density = density./peak;
F = [.05 .1 .05; .1 .4 .1; .05 .1 .05];

for n=1:nsmooth
    
    density = conv2(density,F,'same');
    peak = max(max(density));
    mode = [ X(find(density==peak)), Y(find(density==peak)) ];
    density = density./peak;
    
end

X = [X(1,:); X; X(1,:)];
Y = [ybnd(1)*ones(1,size(Y,2)); Y; ybnd(2)*ones(1,size(Y,2))];
density = [zeros(1,size(density,2)); density; zeros(1,size(density,2))];

X = [xbnd(1)*ones(size(X,1),1), X, xbnd(2)*ones(size(X,1),1)];
Y = [Y(:,1), Y, Y(:,1)];
density = [zeros(size(X,1),1), density, zeros(size(X,1),1)];

density = max(density,0);

nlevels = 20;

levels = 0.01 + (0.7-0.01).*(0:(nlevels-1))./(nlevels-1);
prob = zeros(1,nlevels);
C = contour(X,Y,density,levels);

for i=1:nlevels
    q = find( C(1,:) == levels(i) );
    ncont(i) = length(q);
    
    prob(i) = 0;
    for j=1:ncont(i)
        nvert = C(2,q(j));
        cont{i,j} = C(:,q(j)+(1:nvert));
%         figure(123)
%         plot(cont{i,j}(1,:),cont{i,j}(2,:),'-k','LineWidth',4)
%         hold on
        in = inpoly(data,cont{i,j}');
        prob(i) = prob(i) + sum(sum(in))./length(data);
    end
   %disp([levels(i), prob(i)]);
end


% plot(prob,levels);
% xlabel('Interior Probability')
% ylabel('Density level')


qq = find(prob > 0.2 & prob < 0.9999);
prob = prob(qq);
levels = levels(qq);


level68 = interp1(prob,levels,0.68);
level95 = interp1(prob,levels,0.95);



nlevels = 2;
levels = [level95,level68];
prob = zeros(1,6);

C = contour(X,Y,density,levels);
% disp('--')
for i=1:nlevels
    
    q = find( C(1,:) == levels(i) );
    
    ncont(i) = length(q);
    
    prob(i) = 0;
    for j=1:ncont(i)
        nvert = C(2,q(j));
        cont_out{i,j} = C(:,q(j)+(1:nvert));
        in = inpoly(data,cont_out{i,j}');
        prob(i) = prob(i) + sum(sum(in))./length(data);
    end
%     disp([levels(i), prob(i)]);
end



if abs(prob(1)-0.95)>0.02
    disp('Error95')
end
if abs(prob(2)-0.68)>0.02
    disp('Error68')
end
hold on;
plot(mode(1),mode(2),'+k','MarkerSize',30)
hold off;
contoursout = cont_out;
