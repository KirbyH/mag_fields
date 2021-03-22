clear; 

%% compare speed of sending 1 or several random variables at a time 
nElements = 2e6; 
%{
disp('One by one'); 
tic
Z = zeros(nElements, 1); 
for ii = 1:nElements
    Z(ii) = sampleElement(); 
end
toc
histogram(Z, 'BinWidth', 1); 
%}

disp('Send nElements to sampling function'); 
% much, much faster (~ 20x)
tic
Z = sampleElement(nElements); 
toc

histogram(Z, 'BinWidth', 1, 'normalization', 'pdf'); 
xlabel('PDF [-]'); 
ylabel('Atomic number $Z$ [-]'); 

%% test polynomial fit order
H_dist = importdata('Fe_dist.txt'); 
nFit = 8; 

figure; hold on; 
colororder(jet(nFit)); 
plot(H_dist(:,1), H_dist(:,2), 'x'); 
Legend{1} = 'data'; 

x = linspace(min(H_dist(:,1)), max(H_dist(:,1)), 100); 
for ii = 1:nFit
    a = polyfit(H_dist(:,1), H_dist(:,2), ii); 
    plot(x, polyval(a, x)); 
    Legend{ii+1} = sprintf('Order %i', ii); 
end
legend(Legend); 
% conclusion: 5th order poly seems like enough

%% Plot number flux instead of differential flux
writeEnergyGCR(1);  % show differential flux
load('GCR_data.mat');  % load differential flux information
ylims1 = ylim; 

% number flux conversion
% figure; 
nEl = length(energyGCR(:,1)); 

for ii = 1:nEl
    x = linspace(boundsGCR(ii,1), boundsGCR(ii,2), 100); 
    y = polyval(energyGCR(ii,:), x); 
    x = 10.^x; 
    y = 10.^y .* (x/1e6);  % MULTIPLY BY ENERGY X
    
    yyaxis right; 
    Colors = jet(nEl); 
    loglog(x, y, ':', 'LineWidth', 2, 'Color', Colors(ii,:)); 
    grid on; hold on; 
    xlabel('Energy per proton [eV]'); 
    ylabel('Directional flux [nucleon (m$^2$ s sr)$^{-1}$]'); 
    set(gca, 'YScale', 'log'); 
end
ylims2 = ylim; 
ymax = max([ylims1 ylims2]); 
ymin = min([ylims1 ylims2]); 
ylim([ymin ymax]); 
yyaxis left; 
ylim([ymin ymax]); 

%% Test sampleEnergy.m
load('numFluxGCR.mat');  % we'll use this for plotting
elements = [1 2 6 26]; 
elNames = {'Hydrogen', 'Helium', 'Carbon', 'Iron'}; 

nTests = 1e4; 
for ii = 1:4
    ind = elements(ii); 
    H_eV = sampleEnergy(ones(nTests, 1)*ind);  % sample nTests hydrogen energy vals
    figure; 
    [~, edges] = histcounts(log10(H_eV)); 
    H = histogram(H_eV, 10.^edges, 'Normalization', 'pdf'); 
    set(gca, 'YScale', 'log'); 
    set(gca, 'XScale', 'log');

    % plot GCR Distribution on top of histogram
    hold on; 
    Hmax = max(H.Values); 
    x = linspace(boundsGCR(ind,1), boundsGCR(ind,2), 100); 
    y = polyval(energyGCR(ind,:), x); 
    y = 10.^y; 
    y = y./max(y) * Hmax; 
    loglog(10.^x, y, '--', 'LineWidth', 2); 
    
    xlabel('Energy [eV]'); 
    ylabel('pdf [-]'); 
    title(sprintf('PDF monte-carlo sampling for %s', elNames{ii})); 
end

%% sampling Energy speed test with M
M = 150; 
fprintf('M = %i\n', M); 
tic; 
eV = sampleEnergy(ones(1e5, 1), M); 
toc; 

%% testing powerlog sampling
N = 1e4; 

eV = samplePowerlog(N); 
figure; 
[~, edges] = histcounts(log10(eV)); 
H = histogram(eV, 10.^edges, 'Normalization', 'pdf'); 
set(gca, 'YScale', 'log'); 
set(gca, 'XScale', 'log');

