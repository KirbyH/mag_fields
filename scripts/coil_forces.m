function forces = coil_forces(geom, n_coils, radius)
% calculates the forces from a given geometry on each panel of the array
% and stores the output in a matrix forces [x y z Fx Fy Fz]
% 
% INPUTS (all optional) : 
%     geom : a coil geometry given as an array of coordinate pairs [x y]
%     centered around (0, 0)
%     n_coils : number of coils in the array; default is eight
%     radius : radius of halbach array; default is 5 [m]
% 
% OUTPUTS : 
%     forces : force vector [N] on each panel location (x,y,z). Coordinate
%     locations match with the coil_mp locations returned by
%     create_halbach.m 
% 
% Kirby Heck
% 02/20/2021

% parse variable inputs
if ~exist('geom', 'var')
    geom = importdata('nine_panels.txt');
    geom = geom.data;
end
if ~exist('n_coils', 'var')
    n_coils = 8; 
end
if ~exist('radius', 'var')
    radius = 5;  % default: 5 meters
end

[points, coil_mp, dL] = create_halbach(geom, n_coils, radius); 

end

