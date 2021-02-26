function [coil_forces, coil_mp] = analyze_forces(panel_forces, points)
% Calculates the net force on each hoop (coil) by summing over all the
% panels. Places the force on the center of each hoop, as calculated from
% averaging all of the points. 
% 
% INPUTS : 
%     panel_forces : force vectors on each panel [N x 3]
%     points : [pointsPerCoil x 3 x nCoils] location of end points for each
%       coil (hoop) geometry
% 
% OUTPUTS : 
%     hoop_forces : net force on each coil hoop [n_coils x 3]
%     hoop_mp : midpoint of each coil hoop [n_coils x 3]
% 
% Kirby Heck
% 02/23/2021

dims = size(points); 
n_points = dims(1)-1;  % number of unique points per coil
n_coils = dims(3);  % number of coils

% preallocate
coil_mp = zeros(n_coils, 3); 
coil_forces = zeros(size(coil_mp)); 

for ii = 1:n_coils
    coil_mp(ii,:) = mean(points(1:n_points,:,ii));  % average points, do not include duplicate
    
    start_ind = n_coils*(ii-1) + 1; 
    end_ind = n_coils*ii; 
    % sum forces across all panels on each coil
    coil_forces(ii,:) = sum(panel_forces(start_ind:end_ind, :));  
end

end

