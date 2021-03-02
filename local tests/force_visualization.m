% close all; 
clear; 

geom = coil_geom(1,1,19);
figure(); hold on; 
[points, coil_mp, dL] = create_halbach(geom, 8, 5); 
panel_forces = calc_forces(coil_mp, dL, 1e8, points); 
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
% ellipse
geom = coil_geom(1, 0.1, 11);
figure(); hold on; 
[points, coil_mp, dL] = create_halbach(geom, 8, 5); 
panel_forces = calc_forces(coil_mp, dL, 1e6, points); 
coil_forces = analyze_forces(panel_forces, points, 1);
max(vecnorm(coil_forces,2,2))

