function plot_3DField(B_field, points, SAMPLING)
% Plots a quiver3 of the input magnetic field. If SAMPLING is not equal to
% one, requires uniform grid spacing. 

% f = figure; 
if exist('points', 'var')
    plot_halbach(points); 
end

% plot on every SAMPLING points in x,y,z
if ~exist('SAMPLING', 'var')
    % quiver works best if it only has to handle ~1000 vectors. Make a 
    % guess to the sampling rate: 
    nSide = length(B_field)^(1/3); 
    SAMPLING = floor(nSide/10); 
end

% to sample, reshape, take sampling, then reshape back to a list
Nx = length(unique(B_field(:,1))); 
Ny = length(unique(B_field(:,2))); 
Nz = length(unique(B_field(:,3))); 

if SAMPLING ~= 1
    B_plot = reshape(B_field, [Ny, Nx, Nz, 7]); 
    B_plot = B_plot(1:SAMPLING:end, 1:SAMPLING:end, 1:SAMPLING:end, :); 
    newsize = size(B_plot); 
    newrows = prod(newsize(1:3));  % product of first three dimensions
    B_plot = reshape(B_plot, [newrows, 7]); 
else
    B_plot = B_field; 
end

% normalize vector length
B_plot(:,4:6) = B_plot(:,4:6)./B_plot(:,7); 
B_plot(:,7) = log10(B_plot(:,7));  % take log10 of magnetic field magnitude
minB = -6; maxB = 1; 

% plot with quiver3
q3 = quiver3(B_plot(:,1), B_plot(:,2), B_plot(:,3), ...
    B_plot(:,4), B_plot(:,5), B_plot(:,6)); 
axis equal; 

% zlim([-Zgrid Zgrid]);
% xlim([-Xgrid Xgrid]);
% ylim([-Ygrid Ygrid]);

% COLORS: USE EXTERNAL FUNCTIONS (CREDIT GIVEN IN FCN HEADERS)
colormap(redblue);
c1 = colorbar;
c1.Label.String = '\Delta log_{10}(T)'; 
SetQuiverColor(q3, redblue,'mags',B_plot(:,7),'range',[minB, maxB]);  
caxis([minB, maxB]);
xlabel('$x$ [m]'); ylabel('$y$ [m]'); zlabel('$z$ [m]'); 

end
