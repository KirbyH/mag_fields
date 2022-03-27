function writeNumFlux(plotting)
% Writes a 5th orderpolynomial fit for the number flux of particles
% specified in the setup section with the naming scheme [el]_dist.txt
% 
% INPUTS : 
%     plotting (optional) : plots distributions of number flux
% 
% OUTPUTS : (none - writes a file numFluxGCR.m in cd)
% 
% Kirby Heck 03/13/2021

set(0, 'defaultLegendInterpreter', 'latex'); 
set(0, 'defaultTextInterpreter', 'latex'); 
set(0, 'defaultAxesTickLabelInterpreter', 'latex'); 
set(0, 'defaultLineLineWidth', 1); 

Z = {1, 2, 6, 26}; 
el = {'H', 'He', 'C', 'Fe'}; 

ORDER = 5; 

nEl = length(Z); 
energyGCR = zeros(nEl, ORDER+1); 
boundsGCR = zeros(nEl, 2); 
maxGCR = zeros(nEl, 1); 

for ii = 1:nEl
    ind = Z{ii}; 
    filename = sprintf('%s_dist.txt', el{ii}); 
    dat = importdata(filename); 
    x = dat(:,1)+6;  % data x-axis is in MeV
    y = dat(:,2); 
    % convert from diff flux to numerical flux
    x = 10.^x; 
    y = 10.^y .* x/1e6; % divide by 1e6 because y-section normalized to MeV, x is in eV
    x = log10(x); 
    y = log10(y); 
    
    energyGCR(ind,:) = polyfit(x, y, ORDER); 
    boundsGCR(ind,:) = [min(x), max(x)];  
    % bounds on data provided, +6 is for MeV (FROM DATA)
    
    % calculate maximum value for rejection sampling
    f_prime = polyder(energyGCR(ind,:)); 
    r = roots(f_prime); 
    filter = r>boundsGCR(ind,1) & r<boundsGCR(ind,2);  % filter within domain
    extrema = polyval(energyGCR(ind,:), r(filter)); 
    maxGCR(ind) = 10^(max(extrema));  % stores extremum value of the curve
end

energyGCR(energyGCR==0) = NaN; 
boundsGCR(boundsGCR==0) = NaN; 
maxGCR(maxGCR==0) = NaN; 

% same as writeEnergyGCR.m
if exist('plotting', 'var')
    nElMax = length(energyGCR(:,1)); 
    figure; 
    Colors = jet(nElMax); 
    for ii = 1:nElMax
        x = linspace(boundsGCR(ii,1), boundsGCR(ii,2), 100); 
        loglog(10.^x, 10.^(polyval(energyGCR(ii,:), x)), ...
            'Color', Colors(ii,:)); 
        hold on; grid on; 
        plot(xlim, maxGCR(ii)*[1 1], '--', 'Color', Colors(ii,:)); 
    end
    xlabel('Kinetic Energy [eV]'); 
    ylabel('Differential flux [nucleon (m$^2$ s sr MeV)$^{-1}$]'); 
end

save('numFluxGCR', 'energyGCR', 'boundsGCR', 'maxGCR'); 
end

