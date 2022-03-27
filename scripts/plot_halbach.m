function plot_halbach(points, f1)
% Plots a coil array for an arbitrary point cloud. 
% 
% INPUTS : 
%     points : [pointsPerCoil x 3 x nCoils] 3-dimensional array of coil
%     points 
%     f1 : current figure axes
% 
% OUTPUTS : none
% 
% Kirby Heck 
% 02/18/2021

if exist('f1', 'var')
    figure(f1); 
    hold on;
end

dims = size(points); 
M = dims(1)-1;  % number of panels per coil

if length(dims) == 2 % if plotting one coil only, we need a workaround
    points_3d = zeros([dims 2]); 
    points_3d(:,:,1) = points; 
    points = points_3d;  % little switcheroo if plotting 1 coil
    dims = size(points); 
    nCoils = 1; 
else
    nCoils = dims(3); 
end
pointsX = reshape(points(:,1,:), dims([1 3]));  % select only x-points
pointsY = reshape(points(:,2,:), dims([1 3]));  % likewise for y and z
pointsZ = reshape(points(:,3,:), dims([1 3]));  % [#points x #coils]

arrows = zeros(M*nCoils, 3); 
roots = zeros(M*nCoils, 3); 

for ii = 1:nCoils  % loop through number of coils
    plot3(pointsX(:,ii), pointsY(:,ii), pointsZ(:,ii), ...
        '-k','LineWidth', 2, 'MarkerFaceColor', 'r', 'HandleVisibility', 'off'); 
    hold on; 
    start_ind = (ii-1)*M + 1; 
    end_ind = ii*M; 
    arrows(start_ind:end_ind, :) = diff(points(:,:,ii))/2; % half length
    roots(start_ind:end_ind, :) = points(1:M,:,ii);  % compute vectors for quiver3
end

q = quiver3(roots(:,1), roots(:,2), roots(:,3), ...
    arrows(:,1), arrows(:,2), arrows(:,3)); 
q.HandleVisibility = 'off'; 
q.LineWidth = 2; 
q.Color = 'k'; 
q.AutoScale = 'off';
pause(0.1); % this appears to help
q.NodeChildren(2).Visible = 'off';  % hide quiver stems
view(3); axis equal; grid on; 

% xlabel('x')
% ylabel('y')
% zlabel('z')

end