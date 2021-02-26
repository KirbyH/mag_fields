function [geom] = coil_racetrack(r, w, n_p)
% Creates a racetrack with cosine spacing using parameters r, l, n_p
% 
% INPUTS : 
%     r :  semimajor axis (half of total height)
%     w : semiminor axis (half of total width)
%     n_p : target number of points
%     
% OUTPUTS : 
%     geom : [n_p x 2] matrix of (x,y) coordinate pairs
%     
% Kirby Heck
% 02/25/2021

% the function below breaks if r=w (circle). Workaround: call coil_geom.m
if r==w
    geom = coil_geom(r, w, n_p); 
    return; 
end

% This is going to assemble the racetrack in four pieces with
% ~approximately~ n_p points by leveraging symmetry. 

N = ceil(n_p/8);  % one fourth of points
% initial line, N points
geom(1:N, :) = [w*ones(N,1), linspace(0, r-w, N)']; 
theta = linspace(0, pi/2, N)'; 
geom(N+1:2*N-1, :) = [w*cos(theta(2:end)), w*sin(theta(2:end)) + r-w]; 

geom = [geom; flipud(geom(1:end-1,:)).*[-1 1]];  % flip over y-axis
geom = [geom; flipud(geom(1:end-1,:)).*[1 -1]];  % flip over x-axis

% plot(geom(:,1), geom(:,2)); axis equal; 
end

