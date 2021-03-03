function [B_field] = calc_B(coord, coil_mp, dL, I)
% Kirby Heck and Matt Tuman
% 2/18/21
% This function takes a specified coil configuration and position in
% space and computes the B field at this point
%
% INPUTS:
%   dL: Nx3: Length vector of each panel 
%       dL(N,1) = x direction
%       dL(N,2) = y direction
%       dL(N,3) = z direction
%   coil_mp: Nx3: Midpoints of each panel
%       coil_mp(N, 1) = x position
%       coil_mp(N, 2) = y direction
%       coil_mp(N, 3) = z direction
%   I: 1x1: current in coils (A)
%   coord: 1x3: Coordinates of positions at space where the B field is
%          desired
%       coord(1,1) = x position
%       coord(1,2) = y position
%       coord(1,3) = z position
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
end
if size(coord) ~= [1, 3]
    coord = coord'; 
end
r_mat = coord - coil_mp;  % r-vector, [#midpoints x 3]
r = vecnorm(r_mat,2,2);  % length of each vector [#mp x 1]
r_hat = r_mat./r;
dB = cross(dL, r_hat)./(r.^2);  % contribution due to each mp
B_field = sum(dB) * u_0*I/4/pi;  % summed rows and scaled

end