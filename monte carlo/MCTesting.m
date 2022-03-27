clear; 

%% compare speed of sending 1 or several random variables at a time 
nElements = 2e6; 
%{
disp('One by one'); 
tic
Z = zeros(nElements, 1); 
for ii = 1:nElements
    Z(ii) = sampleElement(); 
end
toc
histogram(Z, 'BinWidth', 1); 
%}

disp('Send nElements to sampling function'); 
% much, much faster (~ 20x)
tic
Z = sampleElement(nElements); 
toc

histogram(Z, 'BinWidth', 1, 'normalization', 'pdf'); 
xlabel('PDF [-]'); 
ylabel('Atomic number $Z$ [-]'); 

%% test polynomial fit order
H_dist = importdata('Fe_dist.txt'); 
nFit = 8; 

figure; hold on; 
colororder(jet(nFit)); 
plot(H_dist(:,1), H_dist(:,2), 'x'); 
Legend{1} = 'data'; 

x = linspace(min(H_dist(:,1)), max(H_dist(:,1)), 100); 
for ii = 1:nFit
    a = polyfit(H_dist(:,1), H_dist(:,2), ii); 
    plot(x, polyval(a, x)); 
    Legend{ii+1} = sprintf('Order %i', ii); 
end
legend(Legend); 
% conclusion: 5th order poly seems like enough

%% Plot number flux instead of differential flux
writeEnergyGCR(1);  % show differential flux
load('GCR_data.mat');  % load differential flux information
ylims1 = ylim; 

% number flux conversion
% figure; 
nEl = length(energyGCR(:,1)); 

for ii = 1:nEl
    x = linspace(boundsGCR(ii,1), boundsGCR(ii,2), 100); 
    y = polyval(energyGCR(ii,:), x); 
    x = 10.^x; 
    y = 10.^y .* (x/1e6);  % MULTIPLY BY ENERGY X
    
    yyaxis right; 
    Colors = jet(nEl); 
    loglog(x, y, ':', 'LineWidth', 2, 'Color', Colors(ii,:)); 
    grid on; hold on; 
    xlabel('Energy per proton [eV]'); 
    ylabel('Directional flux [nucleon (m$^2$ s sr)$^{-1}$]'); 
    set(gca, 'YScale', 'log'); 
end
ylims2 = ylim; 
ymax = max([ylims1 ylims2]); 
ymin = min([ylims1 ylims2]); 
ylim([ymin ymax]); 
yyaxis left; 
ylim([ymin ymax]); 

%% PDF from polynomial fit
load('GCR_data.mat');  % load differential flux information
x = linspace(boundsGCR(1,1), boundsGCR(1,2), 100); 
y = polyval(energyGCR(1,:), x); 
x_log = 10.^x; 
y_log = 10.^y; % .* (x/1e6);  % MULTIPLY BY ENERGY X

figure; 
loglog(x_log, y_log); 
A = trapz(x_log, y_log); 
y_norm = y_log./A; 

figure; 
loglog(x_log, y_norm); 
p = polyfit(x, log10(y_norm), 5); 


%% Test sampleEnergy.m
load('GCR_data.mat');  % we'll use this for plotting
elements = [1 2 6 26]; 
elNames = {'Hydrogen', 'Helium', 'Carbon', 'Iron'}; 

nTests = 1e4; 
for ii = 1:4
    ind = elements(ii); 
    H_eV = sampleEnergy(ones(nTests, 1)*ind);  % sample nTests hydrogen energy vals
    figure; 
    [~, edges] = histcounts(log10(H_eV)); 
    H = histogram(H_eV, 10.^edges, 'Normalization', 'pdf'); 
    set(gca, 'YScale', 'log'); 
    set(gca, 'XScale', 'log');

    % plot GCR Distribution on top of histogram
    hold on; 
    Hmax = max(H.Values); 
    x = linspace(boundsGCR(ind,1), boundsGCR(ind,2), 100); 
    y = polyval(energyGCR(ind,:), x); 
    y = 10.^y; 
    y = y./max(y) * Hmax; 
    loglog(10.^x, y, '--', 'LineWidth', 2); 
    
    xlabel('Energy [eV]'); 
    ylabel('pdf [-]'); 
    title(sprintf('PDF monte-carlo sampling for %s', elNames{ii})); 
end

%% sampling Energy speed test with M
M = 150; 
fprintf('M = %i\n', M); 
tic; 
eV = sampleEnergy(ones(1e5, 1), M); 
toc; 

%% testing powerlog sampling
N = 1e4; 

eV = samplePowerlog(N); 
figure; 
[~, edges] = histcounts(log10(eV)); 
H = histogram(eV, 10.^edges, 'Normalization', 'pdf'); 
set(gca, 'YScale', 'log'); 
set(gca, 'XScale', 'log');

xlabel('Energy [eV]'); 
ylabel('PDF [-]'); 
title('Simple Power-log sampling'); 


%% test rand_shielding.m
r_maj = 2; 
AR = 1; 
r_H = 5; 

geom = coil_racetrack(r_maj, AR, 33); 
[points, coil_mp, dL] = create_halbach(geom, 8, r_H); 

[eff, KE, res] = rand_shielding(points, coil_mp, dL, 'KE', 'powerlog'); %,'seed', 'noseed'); 

% plotting/visualization
figure; 
[~, edges] = histcounts(log10(KE)); 
H = histogram(KE, 10.^edges, 'normalization', 'pdf'); 
set(gca, 'XScale', 'log'); 
set(gca, 'YScale', 'log'); 
title('KE vs PDF of particles'); 

figure; 
hits = KE(logical(res)); 
H1 = histogram(hits, 10.^edges, 'normalization', 'pdf'); 
set(gca, 'XScale', 'log'); 
set(gca, 'YScale', 'log'); 
title('KE vs PDF of hits'); 

figure; 
histogram(KE, 10.^edges); 
hold on; 
histogram(hits, 10.^edges); 
set(gca, 'XScale', 'log'); 
ylabel('Number'); 
xlabel('Energy [eV]'); 

%% test sampleDirVector
r_sphere = 10; 
phi = acos(2*rand(100, 1)-1); %ones(100, 1); 
theta = 2*pi*rand(100, 1); %zeros(100, 1); 
v_hat = sampleDirVector(phi, theta, r_sphere, 1);
v_hat = v_hat*r_sphere*2; 

x = r_sphere.*sin(phi).*cos(theta); 
y = r_sphere.*sin(phi).*sin(theta); 
z = r_sphere.*cos(phi); 

figure; 
quiver3(x, y, z, v_hat(:,1), v_hat(:,2), v_hat(:,3), ...
    'AutoScale', 'off'); 
axis equal; 
bnd = r_sphere; 
xlim([-bnd bnd]); 
ylim([-bnd bnd]); 
zlim([-bnd bnd]); 