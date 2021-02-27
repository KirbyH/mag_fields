function [geom] = coil_racetrack(r, w, n_p, method)
% Creates a racetrack with cosine spacing using parameters r, l, n_p
% 
% INPUTS : 
%     r :  semimajor axis (half of total height)
%     w : semiminor axis (half of total width)
%     n_p : target number of points
%     method : (optional) either 1, 2, or 3. Computational method for 
%       calculating the points on the coil. Suggestions: 
%       1 for general debugging (simplest)
%       2 for calculating forces on one coil (equal spacing panels)
%       3 for calculating magnetic fields (fastest)
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

if ~exist('method', 'var')
    method = 3;  % choose method 3 by default
end

if method == 1
    % METHOD 1: EQUAL NUMBER OF POINTS ON STRAIGHT AND CURVED EDGES
    N = ceil(n_p/8);  % one eighth of points
    % initial line, N points
    geom(1:N, :) = [w*ones(N,1), linspace(0, r-w, N)']; 
    theta = linspace(0, pi/2, N+1)'; 
    geom(N+1:2*N, :) = [w*cos(theta(2:end)), w*sin(theta(2:end)) + r-w]; 
    geom = [geom; flipud(geom(1:end-1,:)).*[-1 1]];  % flip over y-axis
    geom = [geom; flipud(geom(1:end-1,:)).*[1 -1]];  % flip over x-axis

elseif method == 2
    % METHOD 2: APPROXIMATELY EQUALLY SPACED POINTS (RECOMMENDED FOR FORCES)
%     warning('coil_racetrack.m: method 2 selected'); 
    len = r-w; 
    arc = pi*w/2; 
    N = ceil(n_p/4);  % one fourth of points
    N1 = round(N*len/(len+arc));  % proportional number of points to lengths
    N2 = N-N1; 

    geom(1:N1, :) = [w*ones(N1,1), linspace(0, r-w, N1)']; 
    theta = linspace(0, pi/2, N2+1)'; 
    geom(N1+1:N, :) = [w*cos(theta(2:end)), w*sin(theta(2:end)) + r-w]; 
    geom = [geom; flipud(geom(1:end-1,:)).*[-1 1]];  % flip over y-axis
    geom = [geom; flipud(geom(1:end-1,:)).*[1 -1]];  % flip over x-axis

else  % if method is not 1 or 2, then this is the default choice
    % METHOD 3: MINIMAL POINTS (RECOMMENDED FOR FIELDS)
    N = ceil(n_p/2);  % half of number of points
    geom(1,:) = [w, w-r]; 
    geom(2,:) = [w, r-w];  % first two points define line
    theta = linspace(0, pi, N-1)'; 
    geom(3:N, :) = [w*cos(theta(2:end)), w*sin(theta(2:end)) + r-w]; 
    geom = [geom; geom(2:end,:).*[-1 -1]];  % flip over x-axis and y-axis

end


% plot(geom(:,1), geom(:,2)); axis equal; 
end

