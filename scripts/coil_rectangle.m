function [geom] = coil_rectangle(r_maj, AR, n_p)
% Creates a rectangle with aspect ratio AR = h/w, where h and w are the
% half-height and half-width of the rectangle, respectively. The nput
% r_maj is assigned to the height if AR>1, else it is assigned to the
% width. An AR of -1 creates a circle. 
% 
% INPUTS : 
%     r_maj :  semimajor axis (half of total height)
%     AR : Aspect ratio of the racetrack; <1 is wider than it is long
%     n_p : target number of points
%     
% OUTPUTS : 
%     geom : [n_p x 2] matrix of (x,y) coordinate pairs
%     
% Kirby Heck
% 03/24/2021

if AR == -1
    geom = coil_geom(r_maj, -AR, n_p); 
    return;  % skip rest of code
elseif AR > 1
    h = r_maj;  % half height
    w = h/AR;   % half width
else
    w = r_maj;
    h = w*AR; 
end

N = ceil(n_p/2+1); 
N1 = ceil(N*h/(h+w));  % points for sides
N2 = N-N1;  % points for top and bottom
geom(1:N1,:) = [w*ones(N1, 1), linspace(-h, h, N1)']; 
top = linspace(w, -w, N2+1)'; 
geom(N1+1:N, :) = [top(2:end), h*ones(N2, 1)]; 
geom = [geom; -geom(2:end, :)]; 

% figure; 
% plot(geom(:,1), geom(:,2)); axis equal; 
end

