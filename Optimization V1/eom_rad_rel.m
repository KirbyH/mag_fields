function [x_dot] = eom_rad_rel(t, x)
% Equation of motion for ode45 describing the kinematics of a charged
% particle in a magnetic field. 
% 
% INPUTS : 
%     q : charge [C]
%     m : mass [kg]
%     x : initial position and momentum, organized in a [6x1] column vec:  
%       [x; y; z; px; py; z]
% 
% OUTPUTS : 
%     x_dot : differentiated dp/dt and velocity vector, organized in [6x1]
%       column vector: [u; v; w; dpx; dpy; dpz]


p_hat = [x(4), x(5), x(6)];

dL = GL().dL; 
coil_mp = GL().coil_mp; 
I = GL().I; 

B_hat = calc_B(x(1:3)', coil_mp, dL, I); 

gamma = sqrt(norm(p_hat)^2 + 1);

dp = cross(p_hat, B_hat)/gamma;
dr = p_hat/gamma;

x_dot = [dr'; dp'];

end
