% explore relationship between net forces on each coil and the aspect ratio
% for both elliptical and racetrack geometries

Colors = parula(5); 

nTests = 100; 
AR = linspace(1, 8, nTests); 
coilMax = zeros(nTests, 2); 
panMax = zeros(nTests, 2); 
panAvg = zeros(nTests, 2); 
Hradius = 5;  % halbach radius [m]
I = 1;  % assume "non-dimensionalizing"

for ii = 1:nTests
    % racetrack geometry
    geom = coil_racetrack(1, 1/AR(ii), 71); 
    [points, coil_mp, dL] = create_halbach(geom, 8, Hradius); 
    panel_forces = calc_forces(coil_mp, dL, I); 
    coil_forces = analyze_forces(panel_forces, points); 
    norm_panels = panel_forces ./ vecnorm(dL,2,2); 
    panMax(ii,1) = max(vecnorm(norm_panels,2,2)); 
    panAvg(ii,1) = mean(vecnorm(norm_panels,2,2)); 
    coilMax(ii,1) = max(vecnorm(coil_forces,2,2)); 
    
    %ellipse geometry
    geom = coil_geom(1, 1/AR(ii), 71); 
    [points, coil_mp, dL] = create_halbach(geom, 8, Hradius); 
    panel_forces = calc_forces(coil_mp, dL, I); 
    coil_forces = analyze_forces(panel_forces, points); 
    norm_panels = panel_forces ./ vecnorm(dL,2,2); 
    panMax(ii,2) = max(vecnorm(norm_panels,2,2)); 
    panAvg(ii,2) = mean(vecnorm(norm_panels,2,2)); 
    coilMax(ii,2) = max(vecnorm(coil_forces,2,2)); 
end

figure; hold on; grid on; 
colororder(Colors); 
yyaxis left; 
plot(AR, coilMax(:,1)); %, '-', 'Color', Colors(1,:)); 
plot(AR, coilMax(:,2)); %, '-', 'Color', Colors(2,:)); 
ylabel('Maximum net coil force, $F/I^2$ [N A$^{-2}$]'); 

yyaxis right; 
plot(AR, panMax(:,1)); %, '--', 'Color', Colors(1,:)); 
plot(AR, panMax(:,2)); %, '--', 'Color', Colors(2,:)); 
ylabel('Maximum force in any panel, $F/I^2 L$ [N A$^{-2}$ m$^{-1}$]'); 

xlabel('Aspect ratio $h/w$ [-]'); 
legend({'Racetrack', 'Ellipse'}, 'location', 'north'); 
title(sprintf('Forces for array radius = %i m', Hradius)); 