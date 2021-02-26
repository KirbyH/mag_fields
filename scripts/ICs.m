%% Create Sphere For Initial Positions
% Matt Tuman
% 2/25/21

clear all
close all

radius = 10000;
n = 20;
[X,Y,Z] = sphere(n);
X = radius*X;
Y = radius*Y;
Z = radius*Z;

data = [X(:), Y(:), Z(:)];
data = unique(data,'rows');
%% Change into 3d array
for ii = 1:length(data(:,1))
    dir = 0-data(ii,:);
    dir = dir/norm(dir);
    info(ii,:) = [data(ii,:), dir];
end


m = 1.67262e-27;
e =1.6022e-12;
v = sqrt(2*e/m);
v = 1.01e6;
info(:,[4:6]) = info(:,[4:6])*v;


filename = 'nine_panels.txt';
n_coils = 8;
radius = 5;
Xgrid = 20;
Ygrid = 20;
Zgrid = 20;
I = 1e6;
pts_grid = 51;
u_0 = 4*pi*1e-7; % magnetic permeability

% imports
geom = importdata(filename);
geom = geom.data;
coil_points = length(geom(:,1)); 

% ============ intro message ============
disp(['Begin timer: running mag_field.m with ' num2str(pts_grid^3) ' points.']); 
disp(['    X-limits: ' num2str(-Xgrid) ' to ' num2str(Xgrid)]); 
disp(['    Y-limits: ' num2str(-Ygrid) ' to ' num2str(Ygrid)]); 
disp(['    Z-limits: ' num2str(-Zgrid) ' to ' num2str(Zgrid)]); 
tic; 
% ========================================

%% Define X Coordinates for Coil Centers
kk = linspace(0,2*pi,n_coils+1);
cntr = kk(1:end-1); % angles for evenly distributed coils
for ii = 1:n_coils % coordinates for circle centers
    center(:,ii) = [radius*cos(cntr(ii)); radius*sin(cntr(ii)); 0];
end


%% Create Array
cc = 1;
% preallocate 
points = zeros(length(geom), 3, n_coils); 

if mod(n_coils, 2) == 0  %n_coils must be divisible by four
    for ii = 1:n_coils  % ii = coil number
        if cc == 3
            geom = flip(geom);
            cc= 1;
        end
        A = [cos(cntr(ii)) -sin(cntr(ii)) 0;
            sin(cntr(ii)) cos(cntr(ii)) 0;
            0 0 1];
        if mod(ii,2) == 1  % if odd, plot in yz plane
            for jj = 1:coil_points
                points(jj,:,ii) = A*(center(:,1)+[0; geom(jj,1); geom(jj,2)]);
            end
        else  % else plot in xz plane
            for jj = 1:coil_points
                points(jj,:,ii) = A*(center(:,1)+[geom(jj,1); 0; geom(jj,2)]);
            end
        end
        cc = cc+1;
    end
else
    error("Number of coils must be divisible by 4"); 
end

x = linspace(-Xgrid, Xgrid, pts_grid);
y = linspace(-Ygrid, Ygrid, pts_grid);
z = linspace(-Zgrid, Zgrid, pts_grid);


%% Compute coil geometry parameters
% the output here is to start with the points array for the coil geometries
% and end with a dL array (length/direction of each panel on each coil) and
% a midpoint array (midpoint of each panel)
M = coil_points-1; 
dL = zeros(n_coils * M, 3); 
coil_mp = zeros(size(dL)); 

for ii = 1:n_coils  % ii loops through n-coils
    start_ind = (ii-1)*M + 1; 
    end_ind = ii*M; 
    
    dL(start_ind:end_ind, :) = diff(points(:,:,ii)); % first order finite diff
    coil_mp(start_ind:end_ind, :) = points(1:M,:,ii) + dL(start_ind:end_ind,:)/2;
end

GL('I', I);
GL('coil_mp', coil_mp);
GL('dL', dL);



for i = 1:length(info(:,1))
    t_span = [0 .01];
    IC = info(i,:);
    [time, x_dot] = ode45(@eom_rad, t_span, IC);
    record{i} = x_dot(:,[1:3]);
end

figure
for ii = 1:length(record)
    clear x y z
    x = record{ii}(:,1);
    y = record{ii}(:,2);
    z = record{ii}(:,3);
    plot3(x, y, z)
    hold on
end
radius = 5;
plot_halbach(points)
xlim([-radius radius])
ylim([-radius radius])
zlim([-radius radius])
