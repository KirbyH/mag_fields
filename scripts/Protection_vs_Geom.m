%% Optimization Plots with Projectile Simulation
% Matt Tuman and Kirby Heck
% 2/28/2021

%% I WOULD HIGHLY SUGGEST NOT USING 30
% as nTests value. I only did it for high resolution plots

clear all
close all

%% Heat Map With Objective Value the Deflection Rate
    % Halbach Array Radius vs Aspect Ratio of Racetrack
%%{
nTests = 8; 
AR = linspace(1, 10, nTests); 
radius = linspace(3,6, nTests);

for ii = 1:length(AR)
    geom = coil_racetrack(1, 1/AR(ii), 33); 
    for jj = 1:length(radius)
        [points, coil_mp, dL] = create_halbach(geom, 8, radius(jj)); 
        defl_rate(ii,jj) = shielding_rate(points, coil_mp, dL);
    end
end

% Plot the results
figure
p1 = pcolor(AR,radius,defl_rate);
set(p1, 'FaceColor', 'interp'); 
set(p1, 'EdgeColor', 'none'); 
title('Deflection Rates Based on Halbach Radius and Aspect Radius')
xlabel('Radius of Halbach Array [m]')
ylabel('Aspect Ratio (Racetrack) [-]')
colorbar
%}

%% Heat Map With Objective Value the Deflection Rate
    % Halbach Array Radius vs Aspect Ratio of Ellipse
%{
nTests = 8; 
AR = linspace(1, 10, nTests); 
radius = linspace(3,10, nTests);

for ii = 1:length(AR)
    geom = coil_geom(1, 1/AR(ii), 33); 
    for jj = 1:length(radius)
        [points, coil_mp, dL] = create_halbach(geom, 8, radius(jj)); 
        defl_rate(ii,jj) = shielding_rate(points, coil_mp, dL);
    end
end

% Plot the results
figure
pcolor(AR,radius,defl_rate);
title('Deflection Rates Based on Halbach Radius and Aspect Radius')
xlabel('Aspect Ratio (Ellipse) [-]')
ylabel('Radius of Halbach Array [m]')
colorbar
%}


%% Heat Map With Objective Value the Deflection Rate
    % Coil Radius vs Aspect Ratio of Ellipse
%{
nTests = 25; 
AR = linspace(1, 10, nTests); 
radius = linspace(0.25,6, nTests);

for ii = 1:length(AR)
    for jj = 1:length(radius)
        geom = coil_geom(radius(jj), radius(jj)/AR(ii), 21);
        [points, coil_mp, dL] = create_halbach(geom, 8, 5);
        defl_rate(ii,jj) = shielding_rate(points, coil_mp, dL);
    end
end

% Plot the results
figure
pcolor(AR,radius,defl_rate);
title('Deflection Rates Based on Halbach Radius and Aspect Radius')
xlabel('Aspect Ratio (Ellipse) [-]')
ylabel('Major Radius of Coil [m]')
colorbar
%}

%% Heat Map With Objective Value the Deflection Rate
    % Coil Radius vs Aspect Ratio of Racetrack
%%{
nTests = 24; 
AR = linspace(1, 10, nTests); 
radius = linspace(0.25,2, nTests);
FoM = zeros(nTests);  % figure of merit
max_force = zeros(nTests);  
I = 1e7; 

for ii = 1:length(AR)
    for jj = 1:length(radius)
        geom = coil_racetrack(radius(jj), radius(jj)/AR(ii), 21);
        [points, coil_mp, dL] = create_halbach(geom, 8, 5);
        % this line takes a while
        defl_rate(ii,jj) = shielding_rate(points, coil_mp, dL, I);
        
        panel_forces = calc_forces(coil_mp, dL, I); 
        coil_forces = analyze_forces(panel_forces, points); 
        max_force(ii,jj) = max(vecnorm(coil_forces,2,2)); 
        % normalizes to the largest coil force
        FoM(ii,jj) = defl_rate(ii,jj)/max_force(ii,jj); 
    end
end
FoM = FoM./max(FoM, [], 'all'); 

% Plot the results
figure
p3 = pcolor(radius, AR, defl_rate);
set(p3, 'EdgeColor', 'none'); 
set(p3, 'FaceColor', 'interp'); 
title('Deflection Rates Based on Halbach Radius and Aspect Ratio')
xlabel('Major Radius of Coil [m]')
ylabel('Aspect Ratio (Racetrack) [-]')
colormap(redgreen()); 
colorbar
caxis([0 1]); 
%}
    
%% Use fminsearch to do your dirty work for you
%{
for ii = 1:10
    x = rand;
    x(ii,:) = fminsearchbnd(@forFminSearch,[0.5+x*6, 1+x*8.5, 0.1+x*8.5], [0.5, 1, 0.1],[7, 10, 10]);
    display(ii)
end
%}