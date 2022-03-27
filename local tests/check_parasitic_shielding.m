% Calculates the parasitic shielding (amount of radiation that has been
% bent into the spacecraft shield by the magnetic field

set(0, 'defaultLegendInterpreter', 'latex'); 
set(0, 'defaultTextInterpreter', 'latex'); 
set(0, 'defaultAxesTickLabelInterpreter', 'latex'); 
set(0, 'defaultLineLineWidth', 1.5); 

I=4.0e+06; AR=1.50; r_maj=1.85; Hradius=5.6; 
load('rainbow.mat')

geom = coil_racetrack(r_maj, AR, 33); 
% [points, coil_mp, dL] = create_torus(geom, 8, Hradius, 'pumpkin'); 
[points, coil_mp, dL] = create_halbach(geom, 8, Hradius); 

% `sampling` is some sampling "surplus" that corresponds with additional
% area (a solid angle, precisely) to take particles that were "originally
% on a trajectory to avoid the spacecraft" 
sampling = [0:1:2];  % CHANGE THIS TO SAMPLE MORE POINTS (e.g. [0:0.5:10])
nHits = 1000;  % choose N = nHits * sampling so this does not change
nTests = length(sampling); 

nPar = zeros(nTests, 1); 
nDefl = zeros(nTests, 1); 
nInitial = zeros(nTests, 1); 
defl_rate = zeros(nTests, 1); 
par_rate = zeros(nTests, 1); 
ICs = cell(nTests, 1); 
res = cell(nTests, 1); 
res_0 = cell(nTests, 1); 
N = nHits*sampling; 
N(N<nHits) = nHits; 

for ii = 1:nTests
    [ICs{ii}, res{ii}, res_0{ii}] = parasitic_shielding(points, coil_mp, dL, ...
        'N', N(ii), 'I', I, 'sampling', sampling(ii), 'thresh', 1.5, 'randDir', 'thresh'); 
    temp_res = logical(res{ii}); 
    temp_res_0 = logical(res_0{ii}); 
    
    % probably should just do this outside of the for loop
    nPar(ii) = sum(temp_res & ~temp_res_0);  % number of 'parasitic particles'
    nDefl(ii) = sum(~temp_res & temp_res_0);  % number of net deflected particles
    nInitial(ii) = sum(temp_res_0);  % number of particles initially aimed at spacecraft
    defl_rate(ii) = (nDefl(ii)-nPar(ii))/nInitial(ii); 
    par_rate(ii) = nPar(ii)/nDefl(ii); 
end

% plotting
figure; 
colororder(rainbow); 
plot(sampling, defl_rate); 
hold on; 
plot(sampling, par_rate); 
xlabel('Sampling Ratio [-]'); 
ylabel('Performance metric [-]'); 

legend({'True Deflection Rate', 'Parasitic Rate'}, 'location', 'east'); 
grid on; 

%% inward trajectory shielding
% this will take about an hour - you have been warned. 

N = 1e6; 
[IC_in, res_in, res_0_in] = parasitic_shielding(points, coil_mp, dL, ...
    'N', N, 'I', I, 'sampling', sampling(ii), 'thresh', 3, ...
    'randDir', 'inward'); 

%% N-dependence on shielding results for sampling = 4
nHits = 2.^(6:12); 
sampling = 4; 
nTests = length(nHits); 

nPar = zeros(nTests, 1); 
nDefl = zeros(nTests, 1); 
nInitial = zeros(nTests, 1); 
defl_rate = zeros(nTests, 1); 
par_rate = zeros(nTests, 1); 
ICs = cell(nTests, 1); 
res = cell(nTests, 1); 
res_0 = cell(nTests, 1); 
N = nHits*sampling; 

for ii = 1:nTests
    [ICs{ii}, res{ii}, res_0{ii}] = parasitic_shielding(points, coil_mp, dL, ...
        'N', N(ii), 'I', 3e6, 'sampling', sampling, 'thresh', 3, 'randDir', 'thresh'); 
    temp_res = logical(res{ii}); 
    temp_res_0 = logical(res_0{ii}); 
    
    % probably should just do this outside of the for loop
    nPar(ii) = sum(temp_res & ~temp_res_0);  % number of 'parasitic particles'
    nDefl(ii) = sum(~temp_res & temp_res_0);  % number of net deflected particles
    nInitial(ii) = sum(temp_res_0);  % number of particles initially aimed at spacecraft
    defl_rate(ii) = (nDefl(ii)-nPar(ii))/nInitial(ii); 
    par_rate(ii) = nPar(ii)/nDefl(ii); 
end

% plotting
figure; 
plot(N, defl_rate); 
hold on; 
plot(N, par_rate); 

xlabel('Number of samples [-]'); 
ylabel('Effectiveness [-]'); 
title('Convergence of Parasitic Radiation'); 
legend({'Deflection Rate', 'Parasitic Radiation'});


%% Energy Dependence
I=4.0e+06; AR=1.50; r_maj=1.85; Hradius=5.6; 
geom = coil_racetrack(r_maj, AR, 33); 
[points, coil_mp, dL] = create_halbach(geom, 8, Hradius); 

[eV, defl, par] = shielding_vs_energy(r_maj, AR, Hradius, ...
    'N', 3600, 'I', I, 'thresh', 3, 'sampling', 4); 

figure('name', 'Parasitic Radiation vs Energy'); 
colororder(rainbow); 
eV_x = eV*1e-6; 
plot(eV_x, par, eV_x, defl); 
legend({'Parasitic Rate $\xi$','Deflection Rate $\eta_s$'}); %, 'location', 'southoutside', 'box', 'off'); 
set(gca, 'XScale', 'log')
ylim([0 1]); 
xticks([1 3 10 30 100 300 1000]); 
grid on; 
xlabel('Proton Energy [MeV]'); 
ylabel('Performance Metric [-]'); 

%% radial inward only - change in shielding_vs_energy.m
I=4.0e+06; AR=1.50; r_maj=1.85; Hradius=5.6; 
[eV, defl, par] = shielding_vs_energy(r_maj, AR, Hradius, ...
    'I', I, 'N', 3600, 'thresh', 3); 

figure('name', 'Effectiveness vs Energy'); 
colororder(rainbow(2:end,:)); 
eV_x = eV*1e-6; 
plot(eV_x, defl); 
legend('Deflection Rate $\eta_s$'); %, 'location', 'southoutside', 'box', 'off'); 
set(gca, 'XScale', 'log')
ylim([0 1]); 
xticks([1 3 10 30 100 300 1000]); 
grid on; 
xlabel('Proton Energy [MeV]'); 
ylabel('Performance Metric [-]'); 

