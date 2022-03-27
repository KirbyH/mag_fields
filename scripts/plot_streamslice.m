function plot_streamslice(B_field)
% Creates a streamslice plot in the xz (y=0) and xy (z=0) planes as well as
% a 3D plot for the 3D visualization of the magnetic field around the
% Halbach array

Nx = length(unique(B_field(:,1))); 
Ny = length(unique(B_field(:,2))); 
Nz = length(unique(B_field(:,3))); 

X = reshape(B_field(:,1), [Ny, Nx, Nz]); 
Y = reshape(B_field(:,2), [Ny, Nx, Nz]);
Z = reshape(B_field(:,3), [Ny, Nx, Nz]); 
Bx = reshape(B_field(:,4), [Ny, Nx, Nz]); 
By = reshape(B_field(:,5), [Ny, Nx, Nz]); 
Bz = reshape(B_field(:,6), [Ny, Nx, Nz]); 
B = reshape(log10(B_field(:,7)), [Ny, Nx, Nz]);  % log10 of magnetic field

lblue = [0.7 1 1];  % light blue

colorLim = [-4 0]; 
CLab = 'log$_{10} |B|$'; 

% plot slice in y=1 plane
f(1) = figure('name','Field Slices xz'); 
s1 = slice(X, Y, Z, B, [], [1], []); 
axis equal; 
hold on; 
ss1 = streamslice(X,Y,Z,Bx,By,Bz,[],[1],[]); 
set(ss1, 'Color', lblue);
colormap(redblue); 
caxis(colorLim); 
c1 = colorbar; 
c1.Label.String = CLab; 
c1.Label.Interpreter = 'latex'; 
s1.EdgeColor = 'none'; 
s1.FaceColor = 'interp'; 
view(0, 0); 
xlabel('$x$ [m]'); ylabel('$y$ [m]'); zlabel('$z$ [m]'); 
title('Magnetic field visualization for $y=1$'); 

% plot slice in z=0 plane
f(2) = figure('name','Field slices xy'); 
s2 = slice(X, Y, Z, B, [], [], [0]); 
hold on; 
ss2 = streamslice(X,Y,Z,Bx,By,Bz, [], [], [0]);  
set(ss2, 'Color', lblue); 
axis equal; 
colormap(redblue); 
caxis(colorLim); 
c2 = colorbar; 
c2.Label.String = CLab; 
c2.Label.Interpreter = 'latex';
s2.EdgeColor = 'none'; 
s2.FaceColor = 'interp'; 
view([0 90]); 
xlabel('$x$ [m]'); ylabel('$y$ [m]'); zlabel('$z$ [m]'); 
title('Magnetic field visualization for $z=0$'); 

% plot nice 3D figure
f(3) = figure('name','Field 3D streamlines'); hold on; 
s3 = slice(X,Y,Z,B, [max(X(:))], [max(Y(:))], [min(Z(:))]); 
s4 = slice(X,Y,Z,B, [], [], [0]); 
% [Sx, Sy, Sz] = sphere(6);  % start streamlines at radius = 10 [m]
% Sx = Sx*40; Sy = Sy*40; Sz = Sz*40; 
[Sx, Sy, Sz] = meshgrid([-30:20:30]); 
ss3 = streamline(X,Y,Z,Bx,By,Bz,Sx,Sy,Sz); 

set(ss3, 'Color', lblue); 
colormap(redblue); 
caxis(colorLim); 
c3 = colorbar; 
c3.Label.String = CLab; 
c3.Label.Interpreter = 'latex'; 
set(s3, 'EdgeColor', 'none'); 
set(s3, 'FaceColor', 'interp'); 
s4.EdgeColor = 'none'; 
s4.FaceColor = 'interp'; 
view(3); 
xlabel('$x$ [m]'); ylabel('$y$ [m]'); zlabel('$z$ [m]'); 
axis equal; 

end

