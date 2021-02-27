% Checks for convergence of forces depending on the discretization of panel
% geometry

nTests = 30; 
Colors = parula(5); 

nPan = round(linspace(8, 300, nTests));  % number of panels
coilMax = zeros(nTests, 1); 
panMax = zeros(nTests, 1); 
panAvg = zeros(nTests, 1); 
I = 1;  % assume "non-dimensionalizing"

%% TEST PANEL AND COIL FORCES
%%{
tic
for ii = 1:nTests
    % racetrack tests
    geom = coil_racetrack(1, 0.3, nPan(ii)+1, 2); 
    [points, coil_mp, dL] = create_halbach(geom, 8, 5); 
    panel_forces = calc_forces(coil_mp, dL, I); 
    
    coil_forces = analyze_forces(panel_forces, points); 
    norm_panels = panel_forces ./ vecnorm(dL,2,2); 
    panMax(ii,1) = max(vecnorm(norm_panels,2,2)); 
    panAvg(ii,1) = mean(vecnorm(norm_panels,2,2)); 
    coilMax(ii,1) = max(vecnorm(coil_forces,2,2)); 
end

figure; hold on; grid on; 
colororder(Colors); 
yyaxis left 
plot(nPan, coilMax); %, '-', 'Color', Colors(1,:)); 
ylabel('Forces on Coils, $F/I^2$ [N A$^{-2}$]'); 

yyaxis right
plot(nPan, panMax); 
plot(nPan, panAvg); 
ylabel('Forces on Panels, $F/I^2 L$ [N A$^{-2}$ m$^{-1}$]'); 

xlabel('Number of panels $N$ [-]'); 
legend({'Maximum coil force', 'Maximum panel force', 'Averaged panel force'},...
    'location', 'southeast'); 

disp("Computation time for 8-coil Halbach: "); 
toc
%}

%% PANEL FORCES ONLY
% set new values, we can test many more points down here
nTests = 8; 

nPan = round(linspace(10, 10000, nTests));  % number of panels
panMax = zeros(nTests, 1); 
panAvg = zeros(nTests, 1); 

tic
for ii = 1:nTests
    % racetrack tests
    geom = coil_racetrack(1, 0.3, nPan(ii)+1, 2); 
    [points, coil_mp, dL] = create_halbach(geom, 1, 0); 
    panel_forces = calc_forces(coil_mp, dL, I); 
    
    norm_panels = panel_forces ./ vecnorm(dL,2,2); 
    panMax(ii,1) = max(vecnorm(norm_panels,2,2)); 
    panAvg(ii,1) = mean(vecnorm(norm_panels,2,2)); 
end

figure; 
colororder(Colors(2,:)); 
semilogx(nPan, panMax, nPan, panAvg, '--'); 
grid on; 
xlabel('Number of panels $N$ [-]'); 
ylabel('Forces on Panels, $F/I^2 L$ [N A$^{-2}$ m$^{-1}$]'); 
Y = ylim(); 
ylim([0 Y(2)]); 
legend({'Maximum panel force', 'Averaged panel force'}, 'location', 'southeast'); 

disp("Computation time for single panel:"); 
toc
