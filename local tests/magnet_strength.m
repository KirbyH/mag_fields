% Calculates the magnetic field magnitude at the geometric center of
% whatever geometry is given. 

close all; 
Colors = parula(5); 

nTests = 51; 
r = linspace(0.5, 2, nTests);  % semimajor axis, [m] 
AR = linspace(0.8, 5, nTests);  % aspect ratio
I = 1e7; 

r_eff = zeros(nTests, 1); 
AR_eff = zeros(nTests, 1); 

for ii = 1:nTests
    AR_eff(ii) = get_strength(1, AR(ii), I); 
    r_eff(ii) = get_strength(r(ii), 1, I); 
end

figure; grid on; 
colororder(Colors); 
yyaxis left; 
plot(AR_eff, AR); 
ylabel('Aspect Ratio AR [-]'); 

yyaxis right; 
plot(r_eff, AR); 
ylabel('Radius r [m]'); 
xlabel('Magnetic Field Strength $|\vec{B}|$ [T]'); 
xlim([0 50]); 

%% heatmap
%%{
gridB = zeros(nTests, nTests); 
for ii = 1:nTests % loop for each radius
    for jj = 1:nTests  % loop for each AR
        gridB(ii,jj) = get_strength(r(ii), AR(jj), I); 
    end
end
% grid(grid==0) = NaN; 

figure(); 
p = pcolor(AR, r, gridB); hold on; 
[C, h] = contour(AR, r, gridB, [1 5 10 20 50], 'w'); 
set(p, 'EdgeColor', 'none'); 
set(p, 'FaceColor', 'interp'); 
clabel(C, h, 'interpreter', 'latex', 'Color', 'w', 'FontSize', 14); 

hold on; 
xlabel('AR Effects'); 
ylabel('Radius Effects'); 
colormap(hot); 
c1 = colorbar; 
c1.Label.String = 'Strength at center [T]'; 
c1.Label.Interpreter = 'latex'; 
c1.Label.FontSize = 14; 
%}

%% save figures
%%{
savepath = '../figures/'; 
f = findobj('type', 'figure'); 
for k = 1:length(f)
    filename = fullfile(savepath, sprintf('magnet_strength fig %i', k)); 
    saveas(f(k), filename); 
end
%}

%% calculate field strength at center of one coil
function strength = get_strength(r, AR, I)
geom = coil_racetrack(r, r/AR, 73); 
[points, coil_mp, dL] = create_halbach(geom, 1, 0); 
centroid = mean(points);  % [1 x 3] 
strength = vecnorm(calc_B(centroid, coil_mp, dL, I)); 
end

