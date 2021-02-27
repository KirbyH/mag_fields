function [geom] = coil_geom(vert_r, hor_r, n_p)
% Creates geometry for a coil 
% 
% INPUTS : 
%     vert_r : vertical radius of coil [m]
%     hor_r : horizontal radius of coil [m]
%     n_p : number of points to include in coil geometry
% 
% OUTPUTS : 
%     geom : [n_p x 2] coil geometry
%           geom(:,1) = x coordinates
%           geom(:,2) = y coordinates

theta = linspace(0,2*pi, n_p);

for ii = 1:n_p
    r = vert_r*hor_r/sqrt((hor_r*sin(theta(ii)))^2+ (vert_r*cos(theta(ii)))^2);
    x_coord(ii,1) = r*cos(theta(ii));
    y_coord(ii,1) = r*sin(theta(ii));
end

geom = [x_coord,y_coord];

end