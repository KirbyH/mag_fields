function varargout = get_maxForce(points, coil_mp, dL, I)
% Utilizes and combines existing functions calc_forces.m and
% analyze_forces.m to retrieve the maximum net force on any one coil in a
% Halbach array geometry. Basically, I got tired of writing these two
% function calls in scripts. 
% 
% If two output arguments are requested, outputs the maximum compression
% and tension force. Otherwise, just outputs the larger abs. max of the
% two. 
% 
% INPUTS : 
%     usual coil geometry stuff
% 
% OUTPUTS : 
%     maxForce : (1 output argument) absolute maximum coil force
%     maxTension, maxCompression : (2 output arguments) larget forces in
%       both tension and compression. 
% 
% Kirby Heck 
% 03/27/2021

panel_forces = calc_forces(coil_mp, dL, I); 
[net_forces, midpoints] = analyze_forces(panel_forces, points); 

% compressive vs tensile
r_hat = midpoints./vecnorm(midpoints,2,2); 
F_hat = net_forces./vecnorm(net_forces,2,2); 
% determine directions by computing the difference between the coil
% direction and the force direction 
directions = vecnorm(r_hat+F_hat,2,2)-1; 

forces = vecnorm(net_forces,2,2).*directions; 

maxTension = max(forces); 
maxComp = -min(forces(forces<0)); 

if nargout == 1
    varargout = {max([maxTension, maxComp])}; 
else
    varargout = {maxTension, maxComp}; 
end