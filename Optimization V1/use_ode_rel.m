%% Relativistic Script

clear; close all

c = 299792458; % speed of light
m = 1.67262e-27;  % mass of proton [kg]
q = 1.6022e-19;  % charge on a proton, conversion from eV to J 
B_0 = 1;
r_0 = 1;
omega_0 = q*B_0/m;

AR = 0.8688;
C_r = 1.8296;
H_r = 3.9432;

% AR = 2;
% C_r = 0.65;
% H_r = 7;

KE_ev = 1e7;
I = 1e6;

Pos = [1000 -200 0];
vel = -Pos/norm(Pos);

geom = coil_racetrack(C_r/r_0, C_r/AR/r_0, 33); 
[points, coil_mp, dL] = create_halbach(geom, 8, H_r/r_0); 
plots = 1; 


GL('I', I);
GL('coil_mp', coil_mp);
GL('dL', dL);

t = [0 1e-3];
s = omega_0*t;

KE = KE_ev*1.6e-19;
p_hat = sqrt((1+KE/m/c^2)^2-1);
IC  = [Pos p_hat*vel];
[~, traj] = ode45(@eom_rad_rel, s, IC);

% nonRelativistic tracking
m = 1.67262e-27;  % mass of proton [kg]
e = 1.6022e-19;  % charge on a proton, conversion from eV to J 
c = 299792458;  % speed of light, m/s
KE_J = KE_ev*e; 
v = c*sqrt(1-(m*c^2/(KE_J+m*c^2))^2);  % relativistic velocity
IC = [Pos v*vel];
t_span = [0 1e-3];
[~, traj_r] = ode45(@eom_rad, t_span, IC);


plot3(traj(:,1), traj(:,2), traj(:,3))
hold on
plot3(traj_r(:,1), traj_r(:,2), traj_r(:,3))
plot_halbach(points)
legend('Relativistic', 'Non-Relativistic')
xlim([-30 30])
ylim([-30 30])
zlim([-30 30])
title(['Proton Energy: ', num2str(KE_ev/1e6), ' MeV'], ['Current: ', num2str(I/1e6), ' MA'], 'FontSize', 20)
view(3)