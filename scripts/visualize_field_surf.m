function visualize_field_surf(r_maj, AR, r_H, I, options)
% visualizes the magnetic field strength as a function of position in the
% center plane (symmetry plane) of a Halbach array

arguments
    r_maj double = 1.85     % major coil radius for racetrack
    AR double = 1.5         % aspect ratio of racetrack
    r_H double = 5.4        % halbach radius 
    I double = 4e6          % effective coil current
    options.LIM double = r_H*2
    options.zlims double = [-6, 3]
end
   
set(0, 'defaultLegendInterpreter', 'latex'); 
set(0, 'defaultTextInterpreter', 'latex'); 
set(0, 'defaultAxesTickLabelInterpreter', 'latex'); 
set(0, 'defaultColorbarTickLabelInterpreter', 'latex'); 
set(0, 'defaultLineLineWidth', 1.5); 
 
geom = coil_racetrack(r_maj, AR, 51); 
[points, coil_mp, dL] = create_halbach(geom, 8, r_H); %,'pumpkin'); 

xlims = linspace(-options.LIM, options.LIM, 80); 
ylims = xlims; 
zlims = [-6, 3]; 
[X, Y] = meshgrid(xlims, ylims); 
B = zeros(size(X)); 

for ii = 1:length(xlims)
    for jj = 1:length(ylims)
        B(ii,jj) = norm(calc_B([xlims(ii), ylims(jj), 0], coil_mp, dL, I)); 
    end
end
origin_ind = find(X==0&Y==0); 
B = log10(B); 
maxB = max(B, [], 'all'); 
% B = B/maxB;  % this would normalize the field magnitude 
if ~isempty(origin_ind)
    B(origin_ind) = NaN;  % we know this sorta falls apart at (0,0) (fields cancel to zero)
end

% draw cylinder "halo module"
[X_c, Y_c, Z_c] = cylinder(1.5); 
Z_c = 6*(Z_c-0.5);  % stretch Z accordingly

f = figure('name', 'Magnetic Field Visualization'); 
plot_halbach(points, f); 
hold on; 
halo = surf(X_c,Y_c,Z_c); 
colormap('cool')
field = mesh(X,Y,B); 

% make graph pretty
set(halo, 'facecolor', 0.6*[1 1 1]); 
set(halo, 'facealpha', 0.6); 
set(halo, 'edgecolor', 'none'); 
set(field, 'facealpha', 0.5); 
C = colorbar; 
caxis(zlims); 
C.Label.String = 'Field Strength [log(T)]'; 
C.Label.Interpreter = 'latex'; 
zlim(zlims); 

% some lighting options... they kinda suck at the moment
% camlight(45, 45); 
% lighting gouraud

xlabel('$x$ [m]'); 
ylabel('$y$ [m]'); 
zlabel('Field Strength [log$_{10}$(T)]')
end