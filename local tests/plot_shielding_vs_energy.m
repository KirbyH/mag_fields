% Plots a Transmissiblity vs. energy curve for cosmic radiation in a
% Halbach array with `params` configuration. params can be a 2D array to
% compute multiple Halbach parameters. 

clc; clear; %close all; 

% format: [r_maj, AR, Hradius, I]
params = [1.85, 1.5, 5.6, 4e6]; 

%{  
% generates the curve in the final paper
params = [3, -1, 5, 1e6; 
          2, -1, 5, 1e6; 
          1.3, 2, 7, 1e6; 
          1.3, 2, 7, 3e6; 
          1.3, 2, 7, 1e7]; 
%}

Colors = {'r', 'k', 'm', 'b', 'g'}; 
% Colors = strcat(Colors, ':'); 
set(0, 'defaultLineLineWidth', 1.5); 
      
dims = size(params); 
nTests = dims(1); 

figure('name', 'Shielding Curve'); 
hold on; 

for ii = 1:nTests
    r_maj = params(ii,1); 
    AR = params(ii,2); 
    r_H = params(ii,3); 
    I = params(ii,4); 
    [eV, eff(ii,:)] = shielding_vs_energy(r_maj, AR, r_H, ...
        'I', I, 'N', 1000, 'thresh', 3); 
    logEV = log10(eV); 
    semilogx(10.^(logEV-6), 100*(1-eff(ii,:)), [Colors{ii} ':']); 
    hold on; 
    xlabel('Proton energy [MeV]'); 
    ylabel('Transmission probability [\%]'); 
    Legend{ii} = sprintf('$r_H = %.2f, r_{coil} = %.2f, AR = %.2f, I=$%.1E', ...
        r_H, r_maj, AR, I); 
end

legend(Legend, 'location', 'southoutside', 'Box', 'off'); 
set(gca, 'XScale', 'log')
ylim([0 100]); 
xticks([1 3 10 30 100 300 1000]); 
grid on; 

%% replot data with splines
%%{
x = linspace(6, 9, 100); 
x = 10.^x; 

figure; 
for ii = 1:nTests
    y_interp = interp1(eV, 100*(1-eff(ii,:)), x); 
    
    semilogx(x, y_interp, 'Color', Colors{ii}); 
    hold on; 
end
% legend(Legend, 'location', 'southoutside', 'Box', 'off'); 
% set(gca, 'XScale', 'log')
% ylim([0 100]); 
xticks([1 3 10 30 100 300 1000]); 
% grid on; 

%}