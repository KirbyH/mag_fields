%% Current vs Protection Script
clear

%% FMINSEARCH
AR_guess = 1;  
major_r_guess = 0.6; 
Hdist_min_guess = 5; %2;  % this is not the halbach radius, see current_protection.m
I_guess = 0.8e7; 
% thresh = 0;

% for ii = 1:length(thresh)
x{1} = fminsearchbnd(@current_protection,...
    [AR_guess, major_r_guess, Hdist_min_guess, I_guess], ... % initial guesses
    [1, 0.5, 4, 1e5],... % lower bounds
    [5, 2.5, 7,  5e7]); % upper bounds
% end

% fminsearch parameters: (in order)
AR = x{1}(1);
r_maj = x{1}(2);
% hor_r = r_maj*(AR<1) + AR*r_maj*(AR>=1);
% Hradius = x{1}(3) + 1.5 + hor_r;
Hradius = x{1}(3); 
I = x{1}(4);
% thresh = x{1}(5);

%% plot results
I=4.0e+06; AR=1.50; r_maj=1.85; Hradius=5.6; 
% r_maj = 1.35; AR = 1.8; Hradius = 5.85; I = 9.7e6; 
% I=1.1e+07, AR=1.59, r_maj=1.84, Hradius=5.48

geom = coil_racetrack(r_maj, AR, 33);
[points, coil_mp, dL] = create_halbach(geom, 8, Hradius);
figure; 
plot_halbach(points)
defl_rate = rand_shielding(points, coil_mp, dL, ...
    'N', 1024, 'I', I, 'KE', 1e8, 'plots', 'on', 'thresh', 3); 

per = sum(vecnorm(dL,2,2))/8; 
weight = per*(I/400);  % nominal current capacity set at 400 [A]
[maxTens, maxComp] = get_maxForce(points, coil_mp, dL, I); 
maxForce = max([maxTens, maxComp]); 

% defl_rate = parasitic_shielding(points, coil_mp, dL, ...
%     'N', 2048, 'I', I, 'plots', 'on', 'sampling', 4, 'thresh', 1.5); 

fprintf('\n === FINAL RESULTS === \n'); 
fprintf('I=%.1s, AR=%.2f, r_maj=%.2f, Hradius=%.2f\n', I, AR, r_maj, Hradius); 
fprintf('  Deflection rate: %.3f', defl_rate); 
fprintf('  Maximum coil forces - Tension: %.2E N, Compression: %.2E N\n', maxTens, maxComp); 
fprintf('  Length of wire needed: %.2E m\n\n', weight); 


