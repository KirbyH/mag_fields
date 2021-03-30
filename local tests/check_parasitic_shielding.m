% Calculates the parasitic shielding (amount of radiation that has been
% bent into the spacecraft shield by the magnetic field

% r_maj = 2; 
% AR = 1; 
% r_H = 6; 
I=4.0e+06; AR=1.50; r_maj=1.85; Hradius=5.6; 
load('rainbow.mat')

geom = coil_racetrack(r_maj, AR, 33); 
% [points, coil_mp, dL] = create_torus(geom, 8, Hradius, 'pumpkin'); 
[points, coil_mp, dL] = create_halbach(geom, 8, Hradius); 

sampling = 1:0.5:10; 
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

for ii = 1:nTests
    N = nHits*sampling(ii); 
    [ICs{ii}, res{ii}, res_0{ii}] = parasitic_shielding(points, coil_mp, dL, ...
        'N', N, 'I', I, 'sampling', sampling(ii), 'thresh', 3, 'randDir', 'thresh'); 
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
