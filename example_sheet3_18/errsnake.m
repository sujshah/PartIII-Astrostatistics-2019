function [tverts,yverts] = errsnake(tgrid,errorbands)
% outputs the vertices of a polygon representing the error snake at points tgrid of 
% assumes all tgrid is an n-vector, 
% and errorbands is a nx1 or nx2 matrix

%if size(errorbands,2) == 2
    % asymmetric error bands
    tverts = [tgrid;flipud(tgrid)];
    yverts = [errorbands(:,1); flipud(errorbands(:,2))];
%elseif size(errorbands,2) == 1
    % symmetric errorbands
    
%    tverts = [tgrid;flipud(tgrid)];
%    yverts = [errorbands; flipud(errorbands)];
%end



end