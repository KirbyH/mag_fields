function metric = current_protection(x)
%% CURRENTLY SETUP FOR FMINSEARCH
%%% This function computes the ratio of deflection rate to maximum coil force
% INPUTS: 
%           AR: aspect ratio for coil geometry [-]
%           major_r: major radius for coil [m]
%           Hradius_min: minimum halbach array coil distance from any point
%             to the HALO module; 
%             Hradius_min = Hradius - 1.5 [m] - major_r (if AR < 1)
%           I: current in coils [A]
%           type: type of coil (1 = racetrack) (2 = ellipse)
% OUTPUTS:
%           ratio = deflection rate / maximum coil force
% Matt Tuman 3/4/2021

type = 1;
AR = x(1);
r_maj = x(2);
% hor_r = r_maj*(AR<1) + AR*r_maj*(AR>=1);  % calculates the horizontal coil radius
% Hradius = x(3) + hor_r + 1.5;
Hradius = x(3); 
I = x(4);
% thresh = x(5);

% Construct Coil Geometry
if type == 1
    geom = coil_racetrack(r_maj, AR, 33);
else
    geom = coil_geom(r_maj, AR, 21);
end

% Construct array geometry
% [points, coil_mp, dL] = create_torus(geom, 8, Hradius, 'pumpkin');
[points, coil_mp, dL] = create_halbach(geom, 8, Hradius); 

% Compute Deflection Rate
% defl_rate = parasitic_shielding(points, coil_mp, dL, ...
%     'I', I, 'N', 256, 'KE', 1e8, ...'seed', 'noseed', ...
%     'thresh', 1.5, 'sampling', 4); 

defl_rate = rand_shielding(points, coil_mp, dL, ...
    'I', I, 'N', 256, 'KE', 1e8, 'thresh', 3); 

% dims = size(points); 
% nCoils = dims(3); 
per = sum(vecnorm(dL,2,2))/8; 
weight = per*(I/400);  % nominal current capacity set at 400 [A]
[maxTens, maxComp] = get_maxForce(points, coil_mp, dL, I); 
maxForce = max([maxTens, maxComp]); 
% Compute Ratio
alpha = 100/0.4; beta = 100/9e4; lambda = 100/0.5e6; T = 0; 
metric = (defl_rate-T)^-1*alpha + beta*weight + lambda*maxForce; 

fprintf('I=%.1s, AR=%.2f, r_maj=%.2f, Hradius=%.2f\n', I, AR, r_maj, Hradius'); 
fprintf('  Deflection rate: %.3f; metric: %.3f\n', defl_rate, metric); 
fprintf('  Maximum coil forces - Tension: %.2E N, Compression: %.2E N\n', maxTens, maxComp); 
fprintf('  Length of wire needed: %.2E m\n\n', weight); 

end

