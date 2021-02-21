function plot_slices(B_field, points, xslice, yslice, zslice)
% Plot slices of the magnetic field
% 
% INPUTS : 
%     B_field : [N x 7] magnetic field derived from a uniformly spaced grid
%       with columns [ x  y  z  Bx  By  Bz  Bmag ]
%     points (optional) : [pointsPerCoil x 3 x nCoils] plots the Halbach
%     	array if given
%     xslice, yslice, zslice (optional) : points in which to take the
%       visual slices. Default is xslice = [0 10], yslice = [], zslice = [0]
% 
% Kirby Heck
% 02/18/2021

Nx = length(unique(B_field(:,1))); 
Ny = length(unique(B_field(:,2))); 
Nz = length(unique(B_field(:,3))); 

X = reshape(B_field(:,1), [Ny Nx Nz]);  % this Ny Nx ordering is intentional 
Y = reshape(B_field(:,2), [Ny Nx Nz]); 
Z = reshape(B_field(:,3), [Ny Nx Nz]); 
B = reshape(log10(B_field(:,7)), [Ny Nx Nz]); 

% default slice values if not passed arguments
if ~exist('xslice', 'var')
    xslice = [0 10]; 
end
if ~exist('yslice', 'var')
    yslice = []; 
end
if ~exist('zslice', 'var')
    zslice = [0];
end

f = figure; 
% plot Halbach if passed argument
if exist('points', 'var')
    plot_halbach(points, f); 
end

minB = -6; maxB = 1;  % hardcoded colorbar values

s = slice(X, Y, Z, B, xslice, yslice, zslice); 
set(s,'EdgeColor','none');  % hide slice gridlines
set(s, 'FaceColor', 'interp'); 
% set(sl, 'FaceAlpha', 0.4);  % transparency of slices
axis equal; 
colormap(redblue); 
colorbar; 
caxis([minB, maxB]); 

xlabel('x'); ylabel('y'); zlabel('z'); 

end

