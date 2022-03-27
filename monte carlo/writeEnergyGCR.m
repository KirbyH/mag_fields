function writeEnergyGCR(plotting)
% Writes a 5th orderpolynomial fit for the energy distributions in the 
% folder following the naming scheme [el]_dist.txt to a matlab variable 
% file. 
% 
% INPUTS : 
%     plotting (optional) : plots distributions of ion flux
% 
% OUTPUTS : (none - writes a file energyGCR.m in cd)
% 
% Kirby Heck 03/12/2021

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
    dat(:,1) = dat(:,1)+6;  % data x-axis is in MeV
    energyGCR(ind,:) = polyfit(dat(:,1), dat(:,2), ORDER); 
    boundsGCR(ind,:) = [min(dat(:,1)), max(dat(:,1))];  
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

save('GCR_data', 'energyGCR', 'boundsGCR', 'maxGCR'); 
end

