% Compares relativistic to non-relativistic

clear; % close all
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

AR = 1;
C_r = 1.9;
H_r = 5.5;

KE_ev = 1e6; 
I = 1e6; 

% Pos = [1 0 0];
Pos = [100 -50 30]; 
v_hat = -Pos/norm(Pos);

% === generate halbach array coil geometry (dimensional) ===
geom = coil_racetrack(C_r, AR, 33); 
[points, coil_mp, dL] = create_halbach(geom, 8, H_r); 
% ===

% === set global variables for halbach magnetic field calculations ===
GL('I', I);
GL('coil_mp', coil_mp);
GL('dL', dL);
GL('r_0', 1);  % dimensional coordinates, r_0 = 1
% ===

%% non-Relativistic tracking
KE_J = KE_ev*q; 
v = c*sqrt(1-(m*c^2/(KE_J+m*c^2))^2);  % relativistic velocity
% v = sqrt(2*KE_J/m);  % non-relativistic velocity

IC = [Pos, v*v_hat];
t_end = 2*norm(Pos)/norm(v);  % time span to travel across diameter of sphere
t_span = [0 t_end];
[dt, traj_nr] = ode45(@eom_rad, t_span, IC);

%% Relativistic
%%{
% === create halbach coil geometry (non-dimensional) ===
geom = coil_racetrack(C_r/r_0, AR, 33); 
[points, coil_mp, dL] = create_halbach(geom, 8, H_r/r_0); 
% ===

% === update global variables to match ===
GL('I', I);
GL('coil_mp', coil_mp);
GL('dL', dL);
GL('r_0', r_0); 
% ===
%}
GL('r_0', r_0); 

s = omega_0*t_span;  % new timespan scaled by omega_0
p_hat = sqrt((1+KE_J/m/c^2)^2-1);  % KE as given by Tipler and llewellyn
% p_hat = KE_J/m/c^2;  % KE as given by Paolo - incorrect

% opts = odeset('RelTol',1e-4,'AbsTol',1e-6); 
IC  = [Pos/r_0, p_hat*v_hat];  % scaled position ICs
[ds, traj] = ode45(@eom_rad_rel, s, IC);
traj(:, 1:3) = traj(:, 1:3)*r_0;  % rescale trajectory vector

%% plotting - halbach array tests (commented out)
%%{
figure; 
plot3(traj(:,1), traj(:,2), traj(:,3), 'Linewidth', 2)
hold on
plot3(traj_nr(:,1), traj_nr(:,2), traj_nr(:,3), 'rx')

% recreate dimensioned halbach
geom = coil_racetrack(C_r, AR, 33); 
[points, ~, ~] = create_halbach(geom, 8, H_r); 
plot_halbach(points)
legend('Relativistic', 'Non-Relativistic')
bounds = 30; 
xlim(bounds*[-1 1])
ylim(bounds*[-1 1])
zlim(bounds*[-1 1])
title(['Proton Energy: ', num2str(KE_ev/1e6), ' MeV; Current: ', num2str(I/1e6), ' MA']); %, 'FontSize', 20)
view(3)

% plot threshold
%{
thresh = 3; 
[x_s, y_s, z_s] = sphere(10); 
x_s = x_s*thresh; 
y_s = y_s*thresh; 
z_s = z_s*thresh; 
colormap('gray'); 
surf(x_s, y_s, z_s, 'EdgeColor', 'none', 'FaceAlpha', 0.6);  
legend('Relativistic', 'Non-Relativistic', 'Threshold')
%}

%}

%% plotting - simple magnetic field testing
%{
figure; 
plot3(traj(:,1), traj(:,2), traj(:,3)); 
hold on; 
plot3(traj_nr(:,1), traj_nr(:,2), traj_nr(:,3))
hold off; view(0, 90); axis equal; 
legend('Relativistic', 'Non-Relativistic')
title(['Proton Energy: ', num2str(KE_ev/1e6) ' MeV']); 

%}