function eV = samplePowerlog(N, Emin, Emax)
% Samples energy from a simple power-log E^{-1} from E_min to E_max
% 
% INPUTS : 
%     N : Number of samples to take
%     Emin : (optional) minimum energy to sample from [eV], default 1e6 eV
%     Emax : (optional) maximum energy to sample from [eV], default 1e9 eV
%     
% OUTPUTS : 
%     eV : [N x 1] atomic energy, in electron volts, from sampled
%       distributions 
% 
% Kirby Heck
% 03/20/2021

% rng(10);  % set seed

if nargin == 1
    Emin = 1e6; 
    Emax = 1e9; 
end

ksi = rand(N, 1); 

eV = Emin*(Emax/Emin).^ksi; 

end

