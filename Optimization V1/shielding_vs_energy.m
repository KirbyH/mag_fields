function [eV, eff, par] = shielding_vs_energy(r_maj, AR, r_H, varargin)
% Returns values for shielding effectiveness (as determined by
% rand_shielding.m) as a function of energy and returns values. 
% 
% INPUTS : 
%     r_maj : major radius of halbach coils
%     AR : aspect ratio of halbach coils
%     r_H : Halbach radius
%     varargin : same variable input arguments that rand_shielding takes. 
%       SEE: rand_shielding.m
% 
% OUTPUTS : 
%     eV : [N x 1] energies corresponding to rows of effectiveness values
%     eff : [N x 1] effectiveness (1 - %transmission) of the shield
% 
% Kirby Heck 03/19/2021

if nargin < 3
    r_maj = 3; 
    AR = 1; 
    r_H = 5; 
end

geom = coil_racetrack(r_maj, AR, 33); 
[points, coil_mp, dL] = create_halbach(geom, 8, r_H); 

nTests = 13; 
logEV = linspace(6, 9, nTests); 
eV = 10.^logEV;  % log spacing of energy values, discrete bins

% eV = [1 3 10 30 100 300 1000 3000 1e4 3e4 1e5 3e5 1e6] * 1e6; 
% nTests = length(eV); 

eff = zeros(size(eV)); 
par = zeros(size(eV)); 
for ii = 1:nTests
%     eff(ii) = rand_shielding(points, coil_mp, dL, varargin{:}, 'KE', eV(ii)); 
%     eff(ii) = rand_shielding(points, coil_mp, dL, ...
%         'I', I, 'KE', eV(ii), 'N', 1024, 'rel', 'on', 'thresh', 3); 
%     eff(ii) = rel_shielding_rate(points, coil_mp, dL, I, eV(ii)); 
    [eff(ii), par(ii)] = parasitic_shielding(points, coil_mp, dL, ...
        varargin{:}, 'KE', eV(ii)); 
end

%{
% old plotting: 

figure; 
% hold on; 
semilogx(10.^(logEV-6), 100*(1-eff)); 
xlabel('Proton energy [MeV]'); 
ylabel('Transmission probability [\%]'); 
legend(sprintf('$r_H = %.2f, r_{coil} = %.2f, AR = %.2f$', r_H, r_maj, AR), ...
    'location', 'northwest'); 
ylim([0 100]); 
xticks([1 3 10 30 100 300 1000]); 
grid on; 
%}

end