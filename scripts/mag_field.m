%% Create Coil Array Configuration
% Matt Tuman & Kirby Heck
% 2/10/21
% This script calculates the magnetic field around an n-coil Halbach array
% with arbitrary geometry for a given grid resolution. 
% 
% INPUTS: 
%     filename : text file that describes the geometry of the halbach array
%                NOTE: first and last points must match
%     n_coils : number of coils in the Halbach array (multiple of 4)
%     radius : radius of Halbach array
%     Xgrid : span in x-direction to calculate B-field
%     Ygrid : span in y-dir...
%     Zgrid : span in z-dir...
%     I : current in each coil 
%     pts_grid : number of points in each direction for mesh
% 
% WORKING VARIABLES: 
%     ii : iterates through the number of coils
%     xx : iterates through x-grid
%     yy : iterates through y-grid
%     zz : iterates through z-grid
%     points : points located on the Halbach array [points per coil x 3 x n_coils]
%     M : points per coil-1 (for finite difference calculations)
%     coil_mp : midpoints for ALL of the coils [M*n_coils x 3]
%     dL : vector for each "panel" on ALL the coils [M*n_coils x 3]
%     B_field : scaled magnetic field strength at a given point

clear; close all
set(0, 'defaultLegendInterpreter', 'latex'); 
set(0, 'defaultTextInterpreter', 'latex'); 
set(0, 'defaultAxesTickLabelInterpreter', 'latex'); 

filename = 'nine_panels.txt';
n_coils = 8;
radius = 5;
Xgrid = 12;
Ygrid = 12;
Zgrid = 12;
I = 1e6;
pts_grid = 71;
u_0 = 4*pi*1e-7; % magnetic permeability

% imports
geom = importdata(filename);
geom = geom.data;

% ============ intro message ============
disp(['Begin timer: running mag_field.m with ' num2str(pts_grid^3) ' points.']); 
disp(['    X-limits: ' num2str(-Xgrid) ' to ' num2str(Xgrid)]); 
disp(['    Y-limits: ' num2str(-Ygrid) ' to ' num2str(Ygrid)]); 
disp(['    Z-limits: ' num2str(-Zgrid) ' to ' num2str(Zgrid)]); 
tic; 
% ========================================

%% Create Halbach geometry
if mod(n_coils,4)~=0
    warning('Halbach coil should be divisible by four'); 
end
[points, coil_mp, dL] = create_halbach(geom, n_coils, radius); 


%% Compute B Field
x = linspace(-Xgrid, Xgrid, pts_grid);
y = linspace(-Ygrid, Ygrid, pts_grid);
z = linspace(-Zgrid, Zgrid, pts_grid);

B_field = zeros(pts_grid^3, 7); 
for xx = 1:pts_grid
    for yy = 1:pts_grid
        for zz = 1:pts_grid
            coord = [x(xx),y(yy),z(zz)]; 
            r_mat = coord - coil_mp;  % r-vector, [#midpoints x 3]
            r = vecnorm(r_mat,2,2);  % length of each vector [#mp x 1]
            r_hat = r_mat./r; 
            dB = cross(dL, r_hat)./(r.^2);  % contribution due to each mp
            
            ind = yy + pts_grid*(xx-1) + pts_grid^2*(zz-1);  % linear indexing
            % NOTE: this is a weird way to index. First y, then x, then z
            B_field(ind,1:3) = coord; 
            B = sum(dB) * u_0*I/4/pi;  % scaled, summed
            B_field(ind,4:6) = B; 
            B_field(ind,7) = norm(B);  % calculate magnitude for scaling, coloring
        end
    end
end


%% Plot 3d B Field
% plot_3DField(B_field, points, 5); 

%% Plot slices
% plot_slices(B_field, points); 

%% Plot xz-plane and xy-plane streamslice
plot_streamslice(B_field); 

%% Calculate and plot coil forces
coil_forces(coil_mp, dL, I, points); 

%% end timer
toc