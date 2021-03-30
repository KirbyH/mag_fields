function [B_field] = calc_B(coord, coil_mp, dL, I, r_0)
% This function takes a specified coil configuration and position in
% space and computes the B field at this point
%
% INPUTS:
%     coord: [1 x 3] coordinate position to calculate the magnetic field 
%       from; either dimensional or non-dimensionalized to r_0. 
%       coord(1,1) = x position
%       coord(1,2) = y position
%       coord(1,3) = z position
%     dL: Nx3: (optional) Length vector of each panel - DIMENSIONAL
%       dL(N,1) = x direction
%       dL(N,2) = y direction
%       dL(N,3) = z direction
%     coil_mp: (optional) Nx3: Midpoints of each panel - DIMENSIONAL
%       coil_mp(N, 1) = x position
%       coil_mp(N, 2) = y direction
%       coil_mp(N, 3) = z direction
%     I: (optional) 1x1 current in coils (A)
%     r_0 : (optional) Larmor radius scaling factor for non-dimensional
%       coordinates
%
% OUPTUTS:
%   B_field: 1x3: Magnetic Field information
%       B_field(1,1) = x magnetic field component
%       B_field(1,2) = y magnetic field component
%       B_field(1,3) = z magnetic field component
% 
% Kirby Heck and Matt Tuman
% 2/18/21 Updated 03/25/2021

u_0 = 4*pi*1e-7; % magnetic permeability

% import halbach array coil geometry
if nargin == 1
    dL = GL().dL; 
    coil_mp = GL().coil_mp; 
    I = GL().I; 
    r_0 = GL().r_0; 
elseif nargin == 4 && ~exist('r_0', 'var')
    r_0 = 1;  % if neglected, assume r_0 is 1
end

% ensure coordinate is a row vector
if ~isequal(size(coord), [1, 3])
    coord = coord'; 
end
coord = coord*r_0;  % dimensionalize coordinate! 

% compute magnetic field with DIMENSIONAL coordinates
r_mat = coord - coil_mp;  % r-vector, [#midpoints x 3]
r = vecnorm(r_mat,2,2);  % length of each vector [#mp x 1]
r_hat = r_mat./r;
dB = cross(dL, r_hat)./(r.^2);  % contribution due to each mp
B_field = sum(dB) * u_0*I/4/pi;  % summed rows and scaled

end