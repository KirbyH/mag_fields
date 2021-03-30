% takes a closer look at the initial conditions that cause parasitic
% radiaton for a magnetic shield
% 
% Requires three files from parasitic_shielding.m: 
%   ICs [N x 6] 
%   res [N x 1] logical OR numeric tally of hits with shield on
%   res_0 [N x 1] logica OR numeric tally of hits with shield off
% 
% Kirby Heck 
% 03/27/2021

r_maj = 2; 
AR = 1; 
r_H = 6; 

geom = coil_racetrack(r_maj, AR, 33); 
% [points, coil_mp, dL] = create_torus(geom, 8, r_H, 'pumpkin'); 
[points, coil_mp, dL] = create_halbach(geom, 8, r_H); 

% load('inward_r2AR1H6.mat'); 
% ICs = IC_in; 
% res = res_in; 
% res_0 = res_0_in; 
% load('N3600_r2AR1H6_I3e6.mat'); 
% ICs = ICs{10}; 
% res = res{10}; 
% res_0 = res_0{10}; 

KE_eV = 1e8;  % assume something here lol this is not robust
I = 3e6; 
[ICs, res, res_0] = parasitic_shielding(points, coil_mp, dL, ...
    'sampling', 3, 'I', I, 'N', 1024, 'plots', 'on'); 
pause(0.1); 
% defl_rate = rand_shielding(points, coil_mp, dL, ...
%     'I', I, 'N', 1024, 'plots', 'on'); 
pause(0.1); 

% get_maxForce(points, coil_mp, dL, I)

%%

res = logical(res); 
res_0 = logical(res_0); 
filter = res | res_0;  % filter ICs that hit either with shield on or off
IC_subset = ICs(filter, :); 
res_sub = res(filter); 
res_0_sub = res_0(filter); 

%% re-do integrations - setup (copied from parasitic_shielding.m): 
rel = true; 
m = 1.67262e-27;  % mass of proton [kg]
e = 1.6022e-19;  % charge on a proton, conversion from eV to J 
q = e;  % charge of particle [C]
c = 299792458;  % speed of light, m/s
B_0 = 1; 

R = m*c/q/B_0;  % Larmor radius
omega_0 = q*B_0/m; % cyclotron frequency

KE_J = KE_eV*e; 
v = c*sqrt(1-(m*c^2./(KE_J+m*c^2)).^2);  % relativistic velocity to calculate tspan
r_sphere = 50; 

GL('I', I);
GL('coil_mp', coil_mp);
GL('dL', dL);
if rel
    GL('r_0', R); 
else
    GL('r_0', 1); 
end

%% ode45 integrations for trails, if desired
for ii = 1:length(IC_subset) 
    t = [0 2*r_sphere/v]; 
    s = omega_0*t;  
    
    IC_i = IC_subset(ii,:);
    if rel 
        [~, traj] = ode45(@eom_rad_rel, s, IC_i);
        traj = traj*R;  % re-scale position to dimensionalize
    else
        [~, traj] = ode45(@eom_rad, t, IC_i); 
    end
    streaks{ii} = traj(:,1:3);  % xyz trail of points 
end

%% Finish plots
% with this coloring, we will have the following: 
%   1 : successful shielding (green)
%   2 : parasitic radiation (red)
%   3 : unsuccessful shielding (blue)
coloring = res_0_sub + 2*res_sub; 

r_plot = 8;
LineSp = {'k', 'g:', 'r-', 'b--'}; 
LineTh = [0.1, 0.75, 0.5, 0.5];  % set up some linespec options in advance

% === PLOT WITH SHIELD ON === 
f1 = figure();  
for ii = 1:64; %length(IC_subset)  % plot trajectories
    x = streaks{ii}(:,1);
    y = streaks{ii}(:,2);
    z = streaks{ii}(:,3);
    plot3(x, y, z, LineSp{coloring(ii)+1}, 'LineWidth', LineTh(coloring(ii)+1)); 
    hold on;
end

plot_halbach(points, f1)
axis equal 
xlim([-r_plot r_plot])
ylim([-r_plot r_plot])
zlim([-r_plot r_plot])
set(gca, 'Color', [0.5, 0.5, 0.5]); 
view(3)

xlabel('$x$ [m]'); 
ylabel('$y$ [m]');
zlabel('$z$ [m]');

%     Title = sprintf('Shield on for $N$=%i particles, %.3f effective', nRuns, defl_rate);  
%     title(Title); 

%% plot initial conditions
colors = {'ks', 'go', 'rx', 'b^'};  % successful, parasitic, unsuccessful
coloring = res_0_sub + 2*res_sub;  % if not used above

figure; 
for ii = 1:length(colors)
    subs = coloring==(ii-1); 
    x = IC_subset(subs, 1)*R; 
    y = IC_subset(subs, 2)*R; 
    z = IC_subset(subs, 3)*R; 
    plot3(x,y,z, colors{ii}, 'LineWidth', 1, 'MarkerSize', 4); 
    hold on
end
axis equal
grid on
view(3)