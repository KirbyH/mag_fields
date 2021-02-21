% circle parameters
r = 5;
cxyz = [0; 0; 0]; % center
N = randn(3,1); % normal to circle plane
N = N(:)/norm(N);
Q = null(N');
cfun = @(tt) cxyz + r*Q*[cos(tt); sin(tt)];
close all
hold on
xyz = cfun(linspace(0,2*pi,361));
plot3(xyz(1,:),xyz(2,:),xyz(3,:),'b');

% 3D arrows parameters
m = 10; % number of arrows
h = 0.12*r; % height
w = 0.08*r; % width
dir = -1; % 1 anticlock, -1 clock
p = 10; % #circular discretization of the arrow
phi = linspace(0,2*pi,p+1);
Va = [(w/2)*cos(phi);
      (-dir*h)+zeros(size(phi));
      (w/2)*sin(phi)];
Va(:,end+1) = [0; 0; 0];
F = [ 1:p;
      2:p+1;
     (p+2)+zeros(1,p)].';
F2 = 1:p+1; 
Q = [Q, N];
for k=1:m
    tt = 2*pi*k/m;
    R = [cos(tt) -sin(tt)  0;
         sin(tt)  cos(tt)  0;
         0        0        1];
   Vk = cfun(tt)+  Q*R*Va;
   patch('Faces', F, 'Vertices', Vk', 'FaceColor', 'b', 'EdgeColor', 'none');
   patch('Faces', F2, 'Vertices', Vk', 'FaceColor', 'b', 'EdgeColor', 'none');
end
axis equal; 
view(3)