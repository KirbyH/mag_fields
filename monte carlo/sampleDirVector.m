function [v_hat] = sampleDirVector(phi, theta, r_sphere, thresh, rngseed)
% Uses random variables to assign a direction to N initial points given a
% threshold tolerance value and the initial positions
% 
% INPUTS : 
%     phi : [N x 1] altitude angle (measured off the zenith) of each 
%       initial position
%     theta : [N x 1] azimuth angle (measured off +x) of each initial 
%       position
%     r_sphere : radius of sphere of initial positions
%     thresh : threshold radius of the "spaceship" that determines the
%       amount of spread is allowed in each direction vector
%     rngseed : (optional) seed for the random number generation
%     
% OUTPUTS : 
%     v_hat : [N x 3] unit direction vector of the assigned trajectory
% 
% Kirby Heck 
% 03/21/2021

if exist('rngseed', 'var')
    rng(rngseed+1); 
    % pick a consistent seed, but make it different
end

% this script chooses points on a 2D disk and uses the 2D projection of the
% threshold sphere as viewed from the origin to sample velocity directions.
% I am sure there is a better way to do this, but it's what I could come up 
% with. More information at: 
% https://mathworld.wolfram.com/DiskPointPicking.html

psi = asin(thresh/r_sphere);  % maximum deviation angle
rho = thresh*cos(psi);  % maximum deviation radius 
R_prime = r_sphere*cos(psi);  % minimum distance from plate to origin

N = length(phi); 
theta_hat = 2*pi*rand(N,1);  
r_hat = rho*sqrt(rand(N,1)); 

v_hat = [repmat(R_prime, [N,1]), r_hat.*cos(theta_hat), r_hat.*sin(theta_hat)]'; 
for ii = 1:N
    R_Z = [cos(theta(ii)), -sin(theta(ii)), 0; 
           sin(theta(ii)), cos(theta(ii)), 0; 
           0, 0, 1];  
    R_Y = [sin(phi(ii)), 0, -cos(phi(ii)); ...
           0, 1, 0; ...
           cos(phi(ii)), 0, sin(phi(ii))]; 
    v_i = R_Y*v_hat(:,ii); 
    v_i = R_Z*v_i; 
    v_i = -v_i';  % flip and transpose

    v_hat(:,ii) = v_i./vecnorm(v_i,2,2); 
end

v_hat = v_hat';  % final transposition
end

