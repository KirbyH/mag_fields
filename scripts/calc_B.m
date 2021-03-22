function [B_field] = calc_B(coord, coil_mp, dL, I, r_0)
% Kirby Heck and Matt Tuman
% 2/18/21
% This function takes a specified coil configuration and position in
% space and computes the B field at this point
%
% INPUTS:
%   coord: [1 x 3] coordinate position to calculate the magnetic field from
%       coord(1,1) = x position
%       coord(1,2) = y position
%       coord(1,3) = z position
%   dL: Nx3: (optional) Length vector of each panel 
%       dL(N,1) = x direction
%       dL(N,2) = y direction
%       dL(N,3) = z direction
%   coil_mp: (optional) Nx3: Midpoints of each panel
%       coil_mp(N, 1) = x position
%       coil_mp(N, 2) = y direction
%       coil_mp(N, 3) = z direction
%   I: (optional) 1x1: current in coils (A)
%   r_0: (optional) scaling factor for non-dimensionalization
%
% OUPTUTS:
%   B_field: 1x3: Magnetic Field information
%       B_field(1,1) = x magnetic field component
%       B_field(1,2) = y magnetic field component
%       B_field(1,3) = z magnetic field component

u_0 = 4*pi*1e-7; % magnetic permeability

if nargin == 1
    dL = GL().dL; 
    coil_mp = GL().coil_mp; 
    I = GL().I; 
    r_0 = GL().r_0; 
elseif nargin == 4 && ~exist('r_0', 'var')
    r_0 = 1; 
end

%{
% for testing only
B_field = [sin(coord(2))+1, cos(coord(1)), 1]/r_0; 
B_field = [0 0 1]/r_0; 
return; 
%}

if ~isequal(size(coord), [1, 3])
    coord = coord'; 
end
r_mat = coord - coil_mp;  % r-vector, [#midpoints x 3]
r = vecnorm(r_mat,2,2);  % length of each vector [#mp x 1]
r_hat = r_mat./r;
dB = cross(dL, r_hat)./(r.^2);  % contribution due to each mp
B_field = sum(dB) * u_0*I/4/pi / r_0;  % summed rows and scaled

end