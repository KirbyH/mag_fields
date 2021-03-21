% Compares relativistic to non-relativistic

clear; close all
set(0, 'defaultLegendInterpreter', 'latex'); 
set(0, 'defaultTextInterpreter', 'latex'); 
set(0, 'defaultAxesTickLabelInterpreter', 'latex'); 
set(0, 'defaultLineLineWidth', 1); 

%% setup and constants

c = 299792458; % speed of light
m = 1.67262e-27;  % mass of proton [kg]
q = 1.6022e-19;  % charge on a proton, conversion from eV to J 
B_0 = 1;
r_0 = m*c/q/B_0;
omega_0 = q*B_0/m;  % cyclotron frequency

% AR = 0.8688;
% C_r = 1.8296;
% H_r = 3.9432;

AR = 0.92;
C_r = 2.3;
H_r = 4.1;

KE_ev = 1e7; 
I = 1e6; 

Pos = [1 0 0];
% Pos = [100 -50 30]; 
vel = -Pos/norm(Pos);

geom = coil_racetrack(C_r, AR, 33); 
[points, coil_mp, dL] = create_halbach(geom, 8, H_r); 

GL('I', I);
GL('coil_mp', coil_mp);
GL('dL', dL);
GL('r_0', 1); 

%% non-Relativistic tracking
KE_J = KE_ev*q; 
v = c*sqrt(1-(m*c^2/(KE_J+m*c^2))^2);  % relativistic velocity
% v = sqrt(2*KE_J/m); 
IC = [Pos, v*vel];
t_end = 2*norm(Pos)/norm(v);  % guess a time span 
t_span = [0 t_end];
[dt, traj_nr] = ode45(@eom_rad, t_span, IC);

%% Relativistic
%{
geom = coil_racetrack(C_r/r_0, AR, 33); 
[points, coil_mp, dL] = create_halbach(geom, 8, H_r/r_0); 

GL('I', I);
GL('coil_mp', coil_mp);
GL('dL', dL);
GL('r_0', r_0); 
%}

GL('r_0', r_0); 

t = [0 t_end];
s = omega_0*t*2*pi;

p_hat = sqrt((1+KE_J/m/c^2)^2-1);
IC  = [Pos, p_hat*vel];
[ds, traj] = ode45(@eom_rad_rel, s, IC);

%% plotting
%{
plot3(traj(:,1), traj(:,2), traj(:,3))
hold on
plot3(traj_nr(:,1), traj_nr(:,2), traj_nr(:,3), 'rx')

% we keep switching between these
% geom = coil_racetrack(C_r, AR, 33); 
% [points, coil_mp, dL] = create_halbach(geom, 8, H_r); 
plot_halbach(points)
legend('Relativistic', 'Non-Relativistic')
xlim([-30 30])
ylim([-30 30])
zlim([-30 30])
title(['Proton Energy: ', num2str(KE_ev/1e6), ' MeV; Current: ', num2str(I/1e6), ' MA']); %, 'FontSize', 20)
view(3)
%}

%% alt plot
%%{
plot3(traj(:,1), traj(:,2), traj(:,3)); 
hold on; 
plot3(traj_nr(:,1), traj_nr(:,2), traj_nr(:,3))
hold off; view(0, 90); axis equal; 
legend('Relativistic', 'Non-Relativistic')

% circle_theory = 

%}