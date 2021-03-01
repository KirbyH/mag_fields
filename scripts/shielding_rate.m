function defl_rate = shielding_rate(points, coil_mp, dL, I, plots)
% Function version of IC.m which takes a specific point array and a current
% value and calculates the shield effectiveness from a worst-case scenario
% of evenly distributed charged protons. 
% 
% INPUTS : 
%     points : [pointsPerCoil x 3 x n_coils] array prescribing shield
%       geometry
%     coil_mp : midpoints of each panel [nPoints x 3] with each row [x,y,z]
%     dL : [nPoints x 3] vector length and direction of each panel 
%       corresponding with rows in coil_mp 
%     I : (optional) Current value. Default is set to 1e6. 
%     plots : (optional) If included with any value, will plot the
%       trajectories of inbound particles. Takes about 0.75 s extra
% 
% OUTPUTS : 
%     defl_rate : deflection rate of charged particles; this may be changed
%       in the future for more telling metrics by taking random sampling,
%       monte-carlo bombardment, Q-factor, etc.
% 
% INTERNAL PARAMETERS : 
%     r_sphere : starting point of particles
%     n : number of lines of longitude for initial conditions sphere
%     KE : kinetic energy [eV] of radiation particles
%     thresh : threshold radius to check for radiation
%     r_plot : radius to show in plots
% 
% Matt Tuman & Kirby Heck
% 2/25/21

geom = coil_racetrack(1, 1, 33); 
[points, coil_mp, dL] = create_halbach(geom, 8, 6); 
plots = 1; 
if ~exist('I', 'var')
    I = 1e8; 
end

%% Create Sphere For Initial Positions
r_sphere = 100;  % begin particles at 10 km away
n = 15;
[X,Y,Z] = sphere(n);

r_0 = [X(:), Y(:), Z(:)] * r_sphere;  % rearrange and scale sphere
r_0 = unique(r_0,'rows');
nRuns = length(r_0(:,1)); 

% print banner
disp('========= Field Effectiveness =========')
disp(['  Testing field with N=' num2str(nRuns) ' runs']); 
disp('  Timer start'); 
tic; 

%% Calculate velocities
v_hat = -r_0./vecnorm(r_0,2,2);  % initial velocity direction radial inward
KE = 1e8;  % in eV

m = 1.67262e-27;  % mass of proton [kg]
e = 1.6022e-19;  % charge on a proton, conversion from eV to J 
c = 299792458;  % speed of light, m/s
KE_J = KE*e; 

v = c*sqrt(1-(m*c^2/(KE_J+m*c^2))^2);  % relativistic velocity
% v = sqrt(2*KE_J/m);  % non-relativistic velocity
% v = 1e6;
info = [r_0, v_hat*v];

%% Set global vars
GL('I', I);
GL('coil_mp', coil_mp);
GL('dL', dL);

%% Begin ODE45 integrations
res = zeros(nRuns, 1); 
thresh = 7;  % [m]
plotting = false; 

if exist('plots', 'var')
    plotting = true; 
    streaks = cell(nRuns, 1);  % preallocate :)
end

% if plotting, run extra step in if statement
for ii = 1:nRuns
    t_span = [0 2*r_sphere/v];
    IC = info(ii,:);
    [~, traj] = ode45(@eom_rad, t_span, IC);
    trail = traj(:,1:3);  % xyz trail of points
    
    % check to see if the trail intersepts the sphere bounded by 'thresh'
    mags = vecnorm(trail,2,2); 
    [~,ind] = min(mags);  % find index of nearest approach
    if ind ~= length(trail)  % if nearest isn't at the end...
        if mags(ind+1)>mags(ind-1)  % check if next nearest point is ahead 
            res(ii) = does_it_hit(trail(ind,:), trail(ind+1,:), thresh);
        else  % nearest point must be behind
            res(ii) = does_it_hit(trail(ind,:), trail(ind-1,:), thresh);
        end
    else  % nearest point to origin is at the end
        res(ii) = does_it_hit(trail(ind,:), trail(ind-1,:), thresh);
    end
    
    % store trajectories for plotting if true
    if plotting
        streaks{ii} = trail; 
    end
end

defl_rate = (nRuns-sum(res))/nRuns; 

%% Finish plots
if plotting
    figure(); hold on; 
    LineSp = {'g:', 'r-'}; 
    LineTh = [0.5, 0.5];  % set up some linespec options in advance
    
    for ii = 1:nRuns  % plot trajectories
        x = streaks{ii}(:,1);
        y = streaks{ii}(:,2);
        z = streaks{ii}(:,3);
        plot3(x, y, z, LineSp{res(ii)+1}, 'LineWidth', LineTh(res(ii)+1)); 
    end
    
    r_plot = 25;
    plot_halbach(points)
    xlim([-r_plot r_plot])
    ylim([-r_plot r_plot])
    zlim([-r_plot r_plot])
    set(gca, 'Color', [0.5, 0.5, 0.5]); 
    view(3)
    
    xlabel('$x$ [m]'); 
    ylabel('$y$ [m]');
    zlabel('$z$ [m]');
    
    Title = sprintf('Trajectories for $N$=%i particles, %.3f effective', nRuns, defl_rate);  
    title(Title); 
end

toc; 

end