function c = redblue(m)
% REDBLUE    Shades of red and blue color map
% REDBLUE(M), is an M-by-3 matrix that defines a colormap.
% The colors begin with bright blue, range through shades of
% blue to white, and then through shades of red to bright red.
% REDBLUE, by itself, is the same length as the current figure's
% colormap. If no figure exists, MATLAB creates one.
% 
% For example, to reset the colormap of the current figure:
% 
%         colormap(redblue)
% 
% See also HSV, GRAY, HOT, BONE, COPPER, PINK, FLAG, 
% COLORMAP, RGBPLOT.
% 
% Adam Auton, 9th October 2009
% 
% Updated 02/16/2021 by Kirby Heck to extend the range of the reds and
% blues 

if nargin < 1, m = size(get(gcf,'colormap'),1); end

if (mod(m,2) == 0)
    % From [0 0 1] to [1 1 1], then [1 1 1] to [1 0 0];
    m1 = m*0.5;
    r = (0:m1-1)'/max(m1-1,1);
    g = r;
    r = [zeros(m1, 1); r; ones(m1,1); linspace(1, 0.33, m1)'];
    g = [zeros(m1, 1); g; flipud(g); zeros(m1,1)];
    b = flipud(r);
else
    % From [0 0 1] to [1 1 1] to [1 0 0];
    m1 = floor(m*0.5);
    r = (0:m1-1)'/max(m1,1);
    g = r;
    r = [r; ones(m1+1,1)];
    g = [g; 1; flipud(g)];
    b = flipud(r);
end

c = [r g b]; 

