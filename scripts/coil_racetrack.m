function [geom] = coil_racetrack(r, w, n_p)
% Creates a racetrack with cosine spacing using parameters r, l, n_p
% 
% INPUTS : 
%     r :  semimajor axis (half of total height)
%     w : semiminor axis (half of total width)
%     n_p : number of points
%     
% OUTPUTS : 
%     geom : [n_p x 2] matrix of (x,y) coordinate pairs
%     
% Kirby Heck
% 02/25/2021

% This is going to assemble the racetrack in four pieces with
% ~approximately~ n_p points by leveraging symmetry. 

theta_c = atan((r-w)/w); 
N = ceil(n_p/4); 
n_p = N*4-3;  % actual number of points
geom = zeros(n_p, 2);  % preallocate memory

theta = linspace(0, pi/2, N); 
straight_ind = (theta>

end

