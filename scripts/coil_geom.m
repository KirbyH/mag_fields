function [geom] = coil_geom(r_maj, AR, n_p)
% Creates geometry for a coil 
% 
% INPUTS : 
%     r_maj : major radius of coil [m]
%     AR : Aspect ratio (height/width) [-]. AR < 1 is flatter than tall
%     n_p : number of points to include in coil geometry
% 
% OUTPUTS : 
%     geom : [n_p x 2] coil geometry
%           geom(:,1) = x coordinates
%           geom(:,2) = y coordinates
% 
% UPDATED 03/17 to use major radius and aspect ratio instead

if AR > 1
    vert_r = r_maj; 
    hor_r = r_maj/AR; 
else
    hor_r = r_maj; 
    vert_r = r_maj*AR; 
end

theta = linspace(0,2*pi, n_p);

for ii = 1:n_p
    r = vert_r*hor_r/sqrt((hor_r*sin(theta(ii)))^2+ (vert_r*cos(theta(ii)))^2);
    x_coord(ii,1) = r*cos(theta(ii));
    y_coord(ii,1) = r*sin(theta(ii));
end

geom = [x_coord,y_coord];

end