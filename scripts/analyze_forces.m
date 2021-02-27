function [coil_forces, coil_mp] = analyze_forces(panel_forces, points, plots)
% Calculates the net force on each hoop (coil) by summing over all the
% panels. Places the force on the center of each hoop, as calculated from
% averaging all of the points. 
% 
% INPUTS : 
%     panel_forces : force vectors on each panel [N x 3]
%     points : [pointsPerCoil x 3 x nCoils] location of end points for each
%       coil (hoop) geometry
%     plots (optional) : include if wanting to plot forces. Independent of
%       value of plots variable. 
% 
% OUTPUTS : 
%     hoop_forces : net force on each coil hoop [n_coils x 3]
%     hoop_mp : midpoint of each coil hoop [n_coils x 3]
% 
% Kirby Heck
% 02/23/2021

dims = size(points); 
M = dims(1)-1;  % number of unique points per coil
% n_coils = dims(3);  % number of coils

if length(dims) == 2 % if plotting one coil only, we need a workaround
    points_3d = zeros([dims 2]); 
    points_3d(:,:,1) = points; 
    points = points_3d;  % little switcheroo if plotting 1 coil
    dims = size(points); 
    n_coils = 1; 
else
    n_coils = dims(3); 
end

% preallocate
coil_mp = zeros(n_coils, 3); 
coil_forces = zeros(size(coil_mp)); 

for ii = 1:n_coils
    coil_mp(ii,:) = mean(points(1:M,:,ii));  % average points, do not include duplicate
    
    start_ind = M*(ii-1) + 1; 
    end_ind = M*ii; 
    % sum forces across all panels on each coil
    coil_forces(ii,:) = sum(panel_forces(start_ind:end_ind, :));  
end

if exist('plots', 'var')
    hold on; % plot on top of current figure; 
    plot_halbach(points); 
    
    q1 = quiver3(coil_mp(:,1), coil_mp(:,2), coil_mp(:,3), ...
        coil_forces(:,1), coil_forces(:,2), coil_forces(:,3)); 
    set(q1, 'Color', [0.8 0 0]); 
    set(q1, 'LineWidth', 2); 
end

end

