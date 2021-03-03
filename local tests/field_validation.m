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
dataIn = importdata('Coarse_x50.txt', ' ', 9);  % COMSOL data, 50x50x50 case

load('B_fieldx50.mat'); % loads B_field called B_fieldx50
[err50, avg50] = plotDiffs(B_fieldx50, dataIn, 'Coarse Mesh - Same Conditions'); 
load('B_fieldx50err.mat'); % loads B_field called B_fieldx50
[err50, avg50] = plotDiffs(B_fieldx50, dataIn, 'Coarse Mesh - Racetrack height -20\%'); 
load('B_fieldx50err2.mat'); % loads B_field called B_fieldx50
[err50, avg50] = plotDiffs(B_fieldx50, dataIn, 'Coarse Mesh - Circular coils'); 
load('B_fieldx50err3.mat'); % loads B_field called B_fieldx50
[err50, avg50] = plotDiffs(B_fieldx50, dataIn, 'Coarse Mesh - Wrong current magnitude'); 


%% Fine mesh
load('B_fieldx12.mat'); % loads B_field called B_fieldx12
dataIn = importdata('Fine_x12.txt', ' ', 9);  % COMSOL data, same case
[err12, avg12] = plotDiffs(B_fieldx12, dataIn, 'Fine Mesh'); 


%% function definition
function [err, avg] = plotDiffs(B_field, comsol, titlename)
% visualizes and quantifies differences in magnetic field analyses
% B_field [N x 7] from MATLAB
% comsol {1 x 1} structure with [N x 4] data from COMSOL, same sized mesh 
%   and parameters
BIN_W = 0.01;

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
% filter out extremely wrong values along the z-axis of symmetry
filter = db_diff < 8; 

% error quantification
err = sqrt(1/(length(db_comsol)-1)*sum((db_comsol-db_matlab).^2)); 
err = sqrt(1/(length(db_comsol(filter))-1)*sum((db_comsol(filter)-db_matlab(filter)).^2))
avg = mean(db_comsol-db_matlab);  % avg > 0 means COMSOL is high


figure(); 
xmax = ceil(err*4+0.6)/4;  % round xlim to nearest 0.25; add padding
histogram(db_diff(filter), 'Normalization', 'cdf', 'BinWidth', xmax/50); 
% set(gca, 'YScale', 'log'); 
xlabel('Difference in $|\vec{B}|$, log scale [$\log_{10}$(T)]'); 
ylabel('Points within difference threshold, CDF [-]'); 
if exist('titlename', 'var')
    title(sprintf('%s, Error $\\epsilon = $%.3f', titlename, err)); 
else
    title(sprintf('Error $\\epsilon = $%.3f', err)); 
end
xlim([0 xmax]); 
ylims = ylim(); 
% ylim([0.5, ylims(2)]); 
end