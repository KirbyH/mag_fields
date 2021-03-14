close all; 
clear; 

set(0, 'defaultLegendInterpreter', 'latex'); 
set(0, 'defaultTextInterpreter', 'latex'); 
set(0, 'defaultAxesTickLabelInterpreter', 'latex'); 
set(0, 'defaultLineLineWidth', 0.5); 
%%
r = 0.65; 
AR = 2; 
rH = 7; 
nPoints = 71; 
geom = coil_racetrack(r, r/AR,nPoints, 2);
figure(); hold on; 
[points, coil_mp, dL] = create_halbach(geom, 8, rH); 
panel_forces = calc_forces(coil_mp, dL, 1e7, points); 
coil_forces = analyze_forces(panel_forces, points, 1); 
max(vecnorm(coil_forces,2,2))
xlabel('$x$ [m]'); ylabel('$y$ [m]'); zlabel('$z$ [m]'); 


%% next
% racetrack
%{
geom = coil_racetrack(1, 0.5,33);
figure(); hold on; 
[points, coil_mp, dL] = create_halbach(geom, 8, 5); 
panel_forces = calc_forces(coil_mp, dL, 1e6, points); 
coil_forces = analyze_forces(panel_forces, points, 1); 
max(vecnorm(coil_forces,2,2))
%}

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
geom = coil_racetrack(r, r/AR, nPoints, 1);
figure(); hold on; 
[points, coil_mp, dL] = create_halbach(geom, 1, 0); 
panel_forces = calc_forces(coil_mp, dL, 1e7, points); 
coil_forces = analyze_forces(panel_forces, points, 1); 
max(vecnorm(coil_forces,2,2))

xlim([-2 2]); ylim([-2 2]); 
view(-135, 30)
xlabel('$x$ [m]'); ylabel('$y$ [m]'); zlabel('$z$ [m]'); 

%% EFFECTS OF ARRAY ON ONE COIL
geom = coil_racetrack(r, r/AR,nPoints, 1); 
[points, coil_mp, dL] = create_halbach(geom, 8, rH); 
figure; hold on; 
panel_forces = calc_forces(coil_mp, dL, 1e7, points, 'omitPoints'); 
coil_forces = analyze_forces(panel_forces, points); 
max(vecnorm(coil_forces,2,2))
view(-150, 30); 
xlabel('$x$ [m]'); ylabel('$y$ [m]'); zlabel('$z$ [m]'); 

%% Different views of pringle forces
%{
xlim([2, 10]); ylim([2, 10]); 
view(-75, 40); 
filename = fullfile(savepath, sprintf('force_visualization fig %i', 4)); 
% saveas(gca, filename); 
exportgraphics(gca, [filename '.png']); 


xlim([-1 1]); ylim([0.5, 8]); 
view([120, 40]); 
filename = fullfile(savepath, sprintf('force_visualization fig %i', 5)); 
% saveas(gca, filename); 
exportgraphics(gca, [filename '.png']); 

%}


%% save figures
%%{
savepath = '../figures/'; 
f = findobj('type', 'figure'); 
for k = 1:length(f)
    filename = fullfile(savepath, sprintf('force_visualization fig %i', k)); 
%     saveas(f(k), filename); 
    exportgraphics(f(k), [filename '.png']); 
end
%}