% This script is a complete mess of random tests to try  
%   1. get random sphere sampling correct, and 
%   2. get random conical sampling from a random point on a sphere correct
% 
% Too much of my [Kirby's] life has been dedicated to this completely
% incoherent script for any amount of commenting to do it justice, so
% perhaps just let it be. Or run sections of it if you're curious. 
% 
% Additional resources: 
% https://mathworld.wolfram.com/SpherePointPicking.html
% https://mathworld.wolfram.com/DiskPointPicking.html
% 
% Kirby Heck 09/19/2021

%% number of points for a sphere(n)

n = (1:32)';
nPoints = zeros(length(n), 1); 
for ii = 1:length(n)
    [X,Y,Z] = sphere(n(ii));
    r_0 = [X(:), Y(:), Z(:)];  % rearrange sphere
    r_0 = unique(r_0, 'rows'); 
    nPoints(ii) = length(r_0(:,1)); 
end

% sphere(n) will generate n^2-n+2 = (n+1)(n-2) unique points. 


%% random spherical distribution
N = 512; 
theta = pi*(2*rand(N, 1)-1); 
phi = acos((2*rand(N,1)-1)); 

r_sphere = 50; 
z = cos(phi)*r_sphere; 
x = sqrt(r_sphere^2-z.^2).*cos(theta); 
y = sqrt(r_sphere^2-z.^2).*sin(theta); 

figure('name', 'random sphere points'); 
plot3(x,y,z, 'x', 'MarkerSize', 8, 'linewidth', 0.5); 
hold on; 
plot3(0, 0, 0, 'ro', 'LineWidth', 2); axis equal; 


%% extension: randomly generated perturbations
M = 25; 
thresh = 3;  % inner radius
v = 2*r_sphere; 
psi = asin(thresh/r_sphere); 
rho = thresh*cos(psi); 
th_hat = 2*pi*rand(M,1); 
r_hat = rho*sqrt(rand(M,1)); 
R_prime = r_sphere*cos(psi); 
[x_s, y_s, z_s] = sphere(10); 
x_s = x_s*thresh; 
y_s = y_s*thresh; 
z_s = z_s*thresh; 

figure; 
for ii = 1:N
    v_hat = [repmat(R_prime, [M 1]), r_hat.*cos(th_hat), r_hat.*sin(th_hat)]'; 
    R_Z = [cos(theta(ii)), -sin(theta(ii)), 0; 
           sin(theta(ii)), cos(theta(ii)), 0; 
           0, 0, 1]; 
%     alpha = pi/2-phi(ii);  
%     R_Y = [cos(alpha), 0, -sin(alpha); 
%            0, 1, 0; 
%            sin(alpha), 0, cos(alpha)]; 
    R_Y = [sin(phi(ii)), 0, -cos(phi(ii)); ...
           0, 1, 0; ...
           cos(phi(ii)), 0, sin(phi(ii))]; 
%     v_hat = -v_hat; 
    v_hat = R_Y*v_hat; 
    v_hat = R_Z*v_hat; 
    v_hat = -v_hat';  % flip and transpose

    v_i = v*v_hat./vecnorm(v_hat,2,2); 
    plot3(x(ii), y(ii), z(ii), 'bo', 'LineWidth', 2); 
    hold on; 
    plot3(0,0,0, 'rx', 'LineWidth', 2); 
    X = repmat(x(ii), [M 1]); 
    Y = repmat(y(ii), [M 1]); 
    Z = repmat(z(ii), [M 1]); 
    quiver3(X, Y, Z, v_i(:,1), v_i(:,2), v_i(:,3)); 
    colormap(gray); 
    surf(x_s, y_s, z_s, 'EdgeColor', 'none', 'facealpha', 0.7);  
    hold off; axis equal; 
    
    break  % comment this line out if looking at multiple random samples
    
%     god bless this finally works idk why the rotation matrices are
%     flipped in their multiplication order but uh whatever works I guess

%     view(90+theta(ii)*180/pi, 90-phi(ii)*180/pi)
    
end


%% WRONG sphere distribution
N = 10000; 
theta = 2*pi*rand(N, 1); 
phi = 2*pi*rand(N, 1) - pi; 

x = sqrt(1-cos(phi).^2).*cos(theta); 
y = sqrt(1-cos(phi).^2).*sin(theta); 
z = cos(phi); 

figure; 
plot3(x,y,z, 'x', 'MarkerSize', 2); 
