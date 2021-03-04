% script to detect deviations in the magnetic field strength between MATLAB
% integration and COMSOL. 
% 
% Parameters: 
% 8 coil Halbach
% 5 m Halbach radius
% 1.25 m semimajor axis
% 0.5 m semiminor axis 
% Racetrack Geometry
% COMSOL only - dimension are outer dimensions. Thickness is 0.25 m
% 
% First test case - coarse mesh
% X, Y, Z limits [-50, 50] stepsize 1
% 
% Second test case - fine mesh
% X, Y, Z limits [-12, 12] stepsize 0.2

clear; close all; 

% BINS = [0 0.001 0.01 0.1 1 10]; 
%% Coarse mesh
data50 = importdata('Coarse_x50.txt', ' ', 9);  % COMSOL data, 50x50x50 case

%%{
% figure; hold on; 
load('B_fieldx50.mat'); % loads B_field called B_fieldx50
Title{1} = 'Coarse - Same Geometry'; 
[err50, avg50] = plotDiffs(B_fieldx50, data50, Title{1}, 0.5); 

load('B_fieldx50err.mat'); % loads B_field called B_fieldx50
Title{2} = 'Coarse - Racetrack height -20\%'; 
[err, avg50] = plotDiffs(B_fieldx50, data50, Title{2}, 0.5); 

load('B_fieldx50err2.mat'); % loads B_field called B_fieldx50
Title{3} = 'Coarse - Circular coils'; 
[err, avg50] = plotDiffs(B_fieldx50, data50, Title{3}, 0.5); 

load('B_fieldx50err3.mat'); % loads B_field called B_fieldx50
Title{4} = 'Coarse - Wrong current magnitude'; 
[err, avg50] = plotDiffs(B_fieldx50, data50, Title{4}, 0.5); 

load('B_fieldx50avg.mat'); % loads B_field called B_fieldx50
Title{5} = 'Coarse - Average 3D Coil Dimensions'; 
[err, avg50] = plotDiffs(B_fieldx50, data50, Title{5}, 0.5); 
% hold off; 
%}

%% Fine mesh
load('B_fieldx12.mat'); % loads B_field called B_fieldx12
data12 = importdata('Fine_x12.txt', ' ', 9);  % COMSOL data, same case
[err12, avg12] = plotDiffs(B_fieldx12, data12, 'Fine Mesh'); 


%% Plot difference points
load('B_fieldx50.mat'); % loads B_field called B_fieldx50
plotPoints(B_fieldx50, data50, 3*err50, 'Largest $|\vec{B}|$ difference locations, coarse'); 

% fine mesh
plotPoints(B_fieldx12, data12, 3*err12, 'Largest $|\vec{B}|$ difference locations, fine'); 

%% save figures
%{
savepath = '../figures/'; 
f = findobj('type', 'figure'); 
for k = 1:length(f)
    filename = fullfile(savepath, sprintf('field_validation fig %i', k)); 
    saveas(f(k), filename); 
end
%}


%% function definitions
function [err, avg] = plotDiffs(B_field, comsol, titlename, xmax)
% visualizes and quantifies differences in magnetic field analyses
% B_field [N x 7] from MATLAB
% comsol {1 x 1} structure with [N x 4] data from COMSOL, same sized mesh 
%   and parameters

B_comsol = comsol.data; 
B_matlab = B_field(:,[1 2 3 7]); 

% reorder so rows are the same
B_matlab = sortrows(B_matlab, [1 2 3]); 
B_comsol = sortrows(B_comsol, [1 2 3]); 

% take log10
db_comsol = log10(B_comsol(:,4));  
db_matlab = log10(B_matlab(:,4));  

abs_diff = abs(B_comsol(:,4)-B_matlab(:,4)); 
db_diff = abs(db_comsol-db_matlab); 
% percent_diff = abs_diff./B_matlab(:,4);  % this is a hopeless metric
% filter out extremely wrong values along the z-axis of symmetry
filter = db_diff < 8; 

% error quantification
% err = sqrt(1/(length(db_comsol)-1)*sum((db_comsol-db_matlab).^2)); 
err = sqrt(1/(length(db_comsol(filter))-1)*sum((db_comsol(filter)-db_matlab(filter)).^2)); 
avg = mean(db_comsol-db_matlab);  % avg > 0 means COMSOL is high

figure(); 
if ~exist('xmax', 'var')
    xmax = ceil(err*4+1)/4;  % round xlim to nearest 0.25; add padding 0.25
end
histogram(db_diff(filter), 'Normalization', 'cdf', 'BinWidth', xmax/50); 
xlabel('Difference in $|\vec{B}|$, log scale [$\log_{10}$(T)]'); 
ylabel('Points in Difference Threshold, CDF [-]'); 
if exist('titlename', 'var')
    title(sprintf('%s, Error $\\epsilon = $%.3f', titlename, err)); 
else
    title(sprintf('Error $\\epsilon = $%.3f', err)); 
end
xlim([0 xmax]); 
end

function plotPoints(B_field, comsol, thresh, titlename)
% plots color-coded differences, starting with the largest change in
% magnitude

B_comsol = comsol.data; 
B_matlab = B_field(:,[1 2 3 7]); 

% reorder so rows are the same
B_matlab = sortrows(B_matlab, [1 2 3]); 
B_comsol = sortrows(B_comsol, [1 2 3]); 

% take log10
B_matlab(:,4) = log10(B_matlab(:,4));  
B_comsol(:,4) = log10(B_comsol(:,4));  

db_diff = B_matlab; 
db_diff(:,4) = abs(db_diff(:,4)-B_comsol(:,4)); 
db_diff = flipud(sortrows(db_diff, 4)); 

filter1 = db_diff(:,4)<8;  % filter out zeros
db_diff = db_diff(filter1, :);  % take away egregious outliars
if ~exist('thresh', 'var')
    thresh = 0.5;  % decibel threshold to plot
end

figure; 
filter = db_diff(:,4)>thresh; 
db_plot = db_diff(filter,:); 
db_norm = db_plot(:,4)./max(db_plot(:,4)); 
sc = scatter3(db_plot(:,1), db_plot(:,2), db_plot(:,3), ...
    (db_norm*30).^1.2, db_plot(:,4), 'filled');

set(sc, "MarkerFaceAlpha", 0.6)
c1 = colorbar; 
colormap(hot(12));
c1.Label.String = '\Delta log_{10}(T)'; 
caxis([0.5 3.5]); 
view(3); 
xlabel('$x$ [m]'); ylabel('$y$ [m]'); zlabel('$z$ [m]'); 
if exist('titlename', 'var')
    title(titlename); 
end
axis equal; 
set(gca, 'Color', 0.7*[1 1 1]); 

xmax = max(db_diff(:,1)); 
ymax = max(db_diff(:,2)); 
zmax = max(db_diff(:,3)); 
xlim(xmax*[-1 1]); ylim(ymax*[-1 1]); zlim(zmax*[-1 1]); 
end
