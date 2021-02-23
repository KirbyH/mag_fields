function panel_forces = calc_forces(coil_mp, dL, I, points)
% Calculates the forces from a given geometry on each panel of the array
% and stores the output in a matrix forces [x y z Fx Fy Fz]
% 
% INPUTS : 
%     coil_mp : midpoints of each panel [nPoints x 3] with each row [x,y,z]
%     dL : vector length and direction of each panel corresponding with
%       rows in coil_mp%     
%     I : current thru each coil [A]
%     points (optional) : points for the halbach array. If included, this
%       will plot the halbach array and corresponding forces. 
% 
% OUTPUTS : 
%     forces : force vector [N] on each panel location (x,y,z). Coordinate
%     locations match with the coil_mp locations returned by
%     create_halbach.m 
% 
% Kirby Heck
% 02/20/2021

panel_forces = zeros(size(dL)); 
B = zeros(size(dL));  % this is merely for plotting at the end
nPoints = length(panel_forces(:,1)); 

% loop thru each point, omitting the current node
for ii = 1:nPoints
    subset = ones(1, nPoints); 
    subset(ii) = 0;  % omit this point
    subset = logical(subset);  % needs to be a logical array, not double
    mp_subset = coil_mp(subset, :); 
    dL_subset = dL(subset,:); 
    
    B(ii,:) = calc_B(coil_mp(ii,:), mp_subset, dL_subset, I); 
    panel_forces(ii,:) = cross(dL(ii,:), B(ii,:));  % lorentz force F = I LxB
end
panel_forces = panel_forces*I; 

if exist('points', 'var')
    plot_halbach(points); 
    hold on; 
    q3 = quiver3(coil_mp(:,1), coil_mp(:,2), coil_mp(:,3), ...
        panel_forces(:,1), panel_forces(:,2), panel_forces(:,3), 'Color', 'r'); 
    q4 = quiver3(coil_mp(:,1), coil_mp(:,2), coil_mp(:,3), ...
        B(:,1), B(:,2), B(:,3), 'Color', 'b');
end
end

