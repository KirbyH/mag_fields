close all; 
clear; 

set(0, 'defaultLegendInterpreter', 'latex'); 
set(0, 'defaultTextInterpreter', 'latex'); 
set(0, 'defaultAxesTickLabelInterpreter', 'latex'); 
set(0, 'defaultLineLineWidth', 0.5); 
%%
I=4.0e+06; AR=1.50; r_maj=1.85; Hradius=5.6; 
% I = 1e7; AR = 2; r_maj = 0.65; Hradius = 7; 
nPoints = 31; 
geom = coil_racetrack(r_maj, AR,nPoints, 2);
figure('name', 'halbach_forces'); 
hold on; 
[points, coil_mp, dL] = create_halbach(geom, 8, Hradius); 
panel_forces = calc_forces(coil_mp, dL, I, points); 
coil_forces = analyze_forces(panel_forces, points, 1); 
max(vecnorm(coil_forces,2,2))
xlabel('$x$ [m]'); ylabel('$y$ [m]'); zlabel('$z$ [m]'); 
[maxTens, maxComp] = get_maxForce(points, coil_mp, dL, I); 
fprintf('Max Tensile Force: %d, Max compressive force: %d\n', maxTens, maxComp); 

%% one racetrack missing
dims = size(points); 
M = dims(1)-1; 
points_red = points(:,:,1:end-1); 
coil_mp_red = coil_mp(1:end-M,:); 
dL_red = dL(1:end-M,:); 
% points_red = points(:,:,2:end); 
% coil_mp_red = coil_mp(M+1:end,:); 
% dL_red = dL(M+1:end,:); 
figure; 
panel_forces_red = calc_forces(coil_mp_red, dL_red, I, points_red); 
coil_forces_red = analyze_forces(panel_forces_red, points_red, 1); 
set(gcf, 'name', 'Forces for unbalanced array'); 
max(vecnorm(coil_forces_red,2,2))
xlabel('$x$ [m]'); ylabel('$y$ [m]'); zlabel('$z$ [m]'); 

[maxTens, maxComp] = get_maxForce(points_red, coil_mp_red, dL_red, I); 
fprintf('Max Tensile Force: %d, Max compressive force: %d\n', maxTens, maxComp); 

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
geom = coil_racetrack(r_maj, AR, nPoints, 1);
figure(); hold on; 
[points, coil_mp, dL] = create_halbach(geom, 1, 0); 
panel_forces = calc_forces(coil_mp, dL, 1e7, points); 
coil_forces = analyze_forces(panel_forces, points, 1); 
max(vecnorm(coil_forces,2,2))

xlim([-2 2]); ylim([-2 2]); 
view(-135, 30)
xlabel('$x$ [m]'); ylabel('$y$ [m]'); zlabel('$z$ [m]'); 

%% EFFECTS OF ARRAY ON ONE COIL
geom = coil_racetrack(r_maj,AR,nPoints, 1); 
[points, coil_mp, dL] = create_halbach(geom, 8, Hradius); 
figure; hold on; 
panel_forces = calc_forces(coil_mp, dL, 1e7, points, 'omitPoints'); 
coil_forces = analyze_forces(panel_forces, points); 
max(vecnorm(coil_forces,2,2))
view(-150, 30); 
xlabel('$x$ [m]'); ylabel('$y$ [m]'); zlabel('$z$ [m]'); 
% [maxTens, maxComp] = get_maxForce(points, coil_mp, dL, I)

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
%{
savepath = '../figures/'; 
f = findobj('type', 'figure'); 
for k = 1:length(f)
    filename = fullfile(savepath, sprintf('force_visualization fig %i', k)); 
%     saveas(f(k), filename); 
    exportgraphics(f(k), [filename '.png']); 
end
%}