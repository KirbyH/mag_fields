function metric = current_protection(x)
%% CURRENTLY SETUP FOR FMINSEARCH
%%% This function computes the ratio of deflection rate to maximum coil force
% INPUTS: 
%           AR: aspect ratio for coil geometry [-]
%           major_r: major radius for coil [m]
%           Hradius: radius of Halbach array [m]
%           I: current in coils [A]
%           type: type of coil (1 = racetrack) (2 = ellipse)
% OUTPUTS:
%           ratio = deflection rate / maximum coil force
% Matt Tuman 3/4/2021

lambda = 1/100;
type = 0;
AR = x(1);
major_r = x(2);
Hradius = x(3);
I = x(4);
thresh = x(5);

% Compute maximum coil force
% coilMax = force_v_AR_Fun(AR, major_r, Hradius, I, type);

% Construct Coil Geometry
if type == 1
    geom = coil_racetrack(major_r, major_r/AR, 21);
else
    geom = coil_geom(major_r, major_r/AR, 21);
end

% Construct array geometry
[points, coil_mp, dL] = create_halbach(geom, 8, Hradius);

% Compute Deflection Rate
defl_rate = rel_shielding_rate(points, coil_mp, dL, I);
if defl_rate<thresh
    defl_rate = 0;
end
% Compute Ratio
metric = -defl_rate + sum(vecnorm(dL,2,2))*lambda
end

