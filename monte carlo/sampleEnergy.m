function eV = sampleEnergy(Z, M)
% Samples energy from a 5th order polynomial fit for distributions of
% elements from sampleElement.m as given by Simpson 1983. 
% 
% INPUTS : 
%     Z : [N x 1] atomic numbers to spit out energies for
%     Required also: 'GCR_data.mat' in working directory. This file is a
%       matlab workspace containing: 
%       energyGCR [N x 6] - 5th order polynomial fit, in loglog scale
%       boundsGCR [N x 2] - min and max energy for each element to
%       interpolate between, in log scale
%       maxGCR [N x 1] - maximum differential flux for each element, sets
%       upper bound for rejection sampling, in absolute scale
%     M : (optional) integer that takes groupings of random variables
%       instead of one by one to speed up the rejection sampling
%     
% OUTPUTS : 
%     eV : [N x 1] atomic energy, in electron volts, from sampled
%       distributions 
% 
% Kirby Heck
% 03/13/2021

% IF THIS FILE DOES NOT EXIST, run writeEnergyGCR.m
%   will load two variables: energyGCR (polynomial fit) and boundsGCR (max
%   and min energies provided by data; do not extrapolate. 
load('numFluxGCR.mat'); 
rng(10);  % set seed

boundsLN = 10.^boundsGCR;  % linear scale bounds

if ~exist('M', 'var')
    M = 100; 
end

nEl = length(Z); 
eV = zeros(nEl, 1); 
% crude, brute force method 
for ii = 1:nEl
    att = 0; 
    while eV(ii)==0 % rejection sampling
        att = att+1; 
        ksi = rand(M, 1); 
        ind = Z(ii);  % element number is also the index number 
        x = boundsLN(ind,1) + ksi*(range(boundsLN(ind,:)));  % energy guess
        
        eta = rand(M, 1); 
        f = polyval(energyGCR(ind,:), log10(x));  % evaluate function at energy x
        f = 10.^f;  % convert to absolute scale; 
        g = maxGCR(ind); 
        
        % faster search
        success = find(eta < f/g, 1); 
        if ~isempty(success)
            eV(ii) = x(success); 
        end
        
        % old slower search for M=1 only
%         if eta<(f/g)
%             fprintf('Attempts before success: %i \n', att); 
%             eV(ii) = x; 
% %             plot(x, eta, 'go'); 
%         end
    end
end

end

