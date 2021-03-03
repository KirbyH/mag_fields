% Calculates the magnetic field magnitude at the geometric center of
% whatever geometry is given. 

nTests = 51; 
r = linspace(0.5, 3, nTests);  % semimajor axis, [m] 
AR = linspace(0.8, 5, nTests);  % aspect ratio
I = 1e7; 

r_eff = zeros(nTests, 1); 
AR_eff = zeros(nTests, 1); 

for ii = 1:nTests
    AR_eff(ii) = get_strength(1, AR(ii), I); 
    r_eff(ii) = get_strength(r(ii), 1, I); 
end

figure; 
yyaxis left; 
plot(AR_eff, AR); 
ylabel('Aspect Ratio AR [-]'); 

yyaxis right; 
plot(r_eff, AR); 
ylabel('Radius r [m]'); 
xlabel('Magnetic Field Strength $|\vec{B}|$ [T]'); 

%% heatmap
grid = zeros(nTests, nTests); 
for ii = 1:nTests % loop for each radius
    for jj = 1:nTests  % loop for each AR
        grid(ii,jj) = get_strength(r(ii), AR(jj), I); 
    end
end
% grid(grid==0) = NaN; 

figure(); 
p = pcolor(AR, r, grid); 
set(p, 'EdgeColor', 'none'); 
set(p, 'FaceColor', 'interp'); 
hold on; 
xlabel('AR Effects'); 
ylabel('Radius Effects'); 
colorbar

%% calculate field strength at center of one coil
function strength = get_strength(r, AR, I)
geom = coil_racetrack(r, r/AR, 73); 
[points, coil_mp, dL] = create_halbach(geom, 1, 0); 
centroid = mean(points);  % [1 x 3] 
strength = vecnorm(calc_B(centroid, coil_mp, dL, I)); 
end