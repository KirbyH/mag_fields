% close all; 
clear; 

geom = coil_racetrack(1,1/3,51);
figure(); hold on; 
[points, coil_mp, dL] = create_halbach(geom, 8, 4); 
panel_forces = calc_forces(coil_mp, dL, 1e7, points); 
coil_forces = analyze_forces(panel_forces, points, 1); 
max(vecnorm(coil_forces,2,2))

%% next
% racetrack
geom = coil_racetrack(1, 0.5,33);
figure(); hold on; 
[points, coil_mp, dL] = create_halbach(geom, 8, 5); 
panel_forces = calc_forces(coil_mp, dL, 1e6, points); 
coil_forces = analyze_forces(panel_forces, points, 1); 
max(vecnorm(coil_forces,2,2))

%% next
% racetrack
geom = coil_racetrack(1, 0.3, 71, 1);
figure(); hold on; 
[points, coil_mp, dL] = create_halbach(geom, 8, 5); 
panel_forces = calc_forces(coil_mp, dL, 1e7, points); 
coil_forces = analyze_forces(panel_forces, points, 1); 
max(vecnorm(coil_forces,2,2))


%% next
%{
% ellipse
geom = coil_geom(1, 0.3,71);
figure(); hold on; 
[points, coil_mp, dL] = create_halbach(geom, 8, 5); 
panel_forces = calc_forces(coil_mp, dL, 1e6, points); 
coil_forces = analyze_forces(panel_forces, points, 1);
max(vecnorm(coil_forces,2,2))
%}

%% SINGLE COIL
% racetrack
geom = coil_racetrack(1, 0.3, 71, 1);
figure(); hold on; 
[points, coil_mp, dL] = create_halbach(geom, 1, 0); 
panel_forces = calc_forces(coil_mp, dL, 1e7, points); 
coil_forces = analyze_forces(panel_forces, points, 1); 
max(vecnorm(coil_forces,2,2))

xlim([-1 1]); ylim([-1 1]); 

%% EFFECTS OF ARRAY ON ONE COIL
nPoints = 51; 
geom = coil_racetrack(1,0.3,nPoints, 1); 
[points, coil_mp, dL] = create_halbach(geom, 8, 5); 
figure; hold on; 
panel_forces = calc_forces(coil_mp, dL, 1e7, points, 'omitPoints'); 
coil_forces = analyze_forces(panel_forces, points, 'plot'); 
