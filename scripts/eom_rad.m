function [x_dot] = eom_rad(t, x)
% Equation of motion for ode45 describing the kinematics of a charged
% particle in a magnetic field. 
% 
% INPUTS : 
%     q : charge [C]
%     m : mass [kg]
%     x : initial position and velocity, organized in a [6x1] column vec:  
%       [x; y; z; u; v; w]
% 
% OUTPUTS : 
%     x_dot : differentiated accel/velocity vector, organized in [6x1]
%       column vector: [u; v; w; ax; ay; az]

q = 1.602e-19;  % charge of a proton
m = 1.67e-27;  % mass of a proton

% input global values from global function
dL = GL().dL; 
coil_mp = GL().coil_mp; 
I = GL().I; 

B_field = calc_B(x(1:3)', coil_mp, dL, I); 
Vel = x(4:6);
a = q/m*cross(Vel, B_field);
x_dot = [Vel; a'];

end
