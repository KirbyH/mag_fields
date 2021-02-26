function hit = does_it_hit(p1, p2, radius)
% This function tells you if a line defined by p1 and p2 goes through a 
% spehre centered at the origin with a specified radius
% Input: p1 = 3x1 or 1x3 and specifies the first point of the line
%        p2 = 3x1 or 1x3 and specifies the second point of the line
%        radius = scalar and specifies the radius of the sphere centered at
%                 the origin
% Output: hit = 1 if the line intersects the sphere
%         hit = 0 if the line does not intersect the sphere


center = 0;
r = radius;

a = norm(p2-p1);
b = 2*((p2(1)-p1(1))*p1(1) + (p2(2)-p1(2))*p1(2) + (p2(3)-p1(3))*(p1(3)));
c = p1(1)^2 + p1(2)^2 + p1(3)^2 - r^2;

solution = b*b - 4*a*c;
if solution<0
    hit = 0;
else
    hit = 1;
end

end
