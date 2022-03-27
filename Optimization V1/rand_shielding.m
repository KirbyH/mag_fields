function [defl_rate, KE, res] = rand_shielding(points, coil_mp, dL, varargin)
% Function version of IC.m which takes a specific point array and a current
% value and calculates the shield effectiveness from a worst-case scenario
% of evenly distributed charged protons. 
% 
% USEAGE : 
%     rand_shielding(points, coil_mp, dL) calculates the shielding
%       rate of the specified wire geometry with N=1000 particles for
%       current I=1e6 A, KE=1 MeV, relativistic EOM, randomly sampled
%       initial positions, radially inward initial velocities, threshold
%       value of 3 [m]
%     rand_shielding(points, coil_mp, dL, ...) allows the specification of
%       any of the parameters below (see VARIABLE INPUTS)
% 
% INPUTS : 
%     Variable input arguments after points, coil_mp, and dL: 
%     points : [pointsPerCoil x 3 x n_coils] array prescribing shield
%       geometry
%     coil_mp : midpoints of each panel [nPoints x 3] with each row [x,y,z]
%     dL : [nPoints x 3] vector length and direction of each panel 
%       corresponding with rows in coil_mp 
% 
% VARIABLE INPUTS : 
%     I : current [A], default is 1e6 [A]
%     KE : Kinetic energy in eV, default is 1e8 eV. Use 'powerlog' for 
%       power-log distrubtion
%     plots : controls plotting, default is 'off'
%     rel : relativistic EOM, default is 'on'
%     randPos : randomly generated initial conditions, default is 'on'
%     randDir : randomly selects a trajectory that will hit the spacecraft
%       if left undeflected, default is 'off'
%     N : number of particles, default is 1000
%     thresh : size of spaceship, default is 1.5 [m]
%     seed : set seed, default is 10. Use 'noseed' for random seed. 
%     
% OUTPUTS : 
%     defl_rate : deflection rate divided by total number of particles
%       (which MAY NOT be equal to N)
%     KE_eV : energies at which particles were shot at
%     res : hit/miss for each respective energy
% 
% INTERNAL PARAMETERS : 
%     r_sphere : starting point of particles
%     r_plot : radius to show in plots
% 
% Matt Tuman & Kirby Heck
% 3/20/21

%% PARSE VARIABLES
if nargin < 3
    error('Not enough input arguments')
end

p = inputParser; 
onoff = {'on','off'};
checkonoff = @(x) any(validatestring(x,onoff));
checkvalidKE = @(x) isnumeric(x) || isequal(x, 'powerlog'); 
checkvalidseed = @(x) isnumeric(x) || isequal(x, 'noseed'); 

addParameter(p, 'I', 1e6, @isnumeric)
addParameter(p, 'KE', 1e8, checkvalidKE)
addParameter(p, 'plots', 'off', checkonoff)
addParameter(p, 'rel', 'on', checkonoff)
addParameter(p, 'randPos', 'on', checkonoff)
addParameter(p, 'randDir', 'off', checkonoff)
addParameter(p, 'N', 1000, @isnumeric)
addParameter(p, 'thresh', 1.5, @isnumeric)
addParameter(p, 'seed', 10, checkvalidseed)

parse(p, varargin{:}); 

% set variables
truefalse = @(onoff) isequal(onoff, 'on'); 
I = p.Results.I; 
KE = p.Results.KE; 
plots = truefalse(p.Results.plots); 
rel = truefalse(p.Results.rel); 
randPos = truefalse(p.Results.randPos); 
randDir = truefalse(p.Results.randDir); 
N = p.Results.N; 
thresh = p.Results.thresh; 

% set seed if requested; default is 10
if isnumeric(p.Results.seed)
    rng(p.Results.seed); 
end

%% testing setup (uncomment)
% AR = 0.8688;
% C_r = 1.8296;
% H_r = 3.9432;
% 
% geom = coil_racetrack(C_r, AR, 33); 
% [points, coil_mp, dL] = create_halbach(geom, 8, H_r); 

%% Create Sphere For Initial Positions
r_sphere = 50;  % begin particles at 50 m away

if randPos
    theta = 2*pi*rand(N, 1); 
    phi = acos(2*rand(N,1)-1); 
    
    r_0(:,1) = sin(phi).*cos(theta); 
    r_0(:,2) = sin(phi).*sin(theta); 
    r_0(:,3) = cos(phi); 
else % evenly spaced points on sphere
    nSphere = ceil(sqrt(N));
    [X,Y,Z] = sphere(nSphere);

    r_0 = [X(:), Y(:), Z(:)];  % rearrange sphere
    r_0 = unique(r_0,'rows');
    % r_0 = r_0(randperm(size(r_0, 1)), :);
    theta = atan2(r_0(:,2), r_0(:,1));  % theta values 
    phi = real(acos(r_0(:,3)./vecnorm(r_0,2,2)));  % phi values; need these for random directions
end
r_0 = r_0*r_sphere; 
nRuns = length(r_0(:,1));  % THIS MAY BE DIFFERENT FROM N

% print banner
disp('========= Field Effectiveness =========')
disp(['  Testing field with N=' num2str(nRuns) ' runs']); 
disp('  Timer start'); 
tic; 

%% Calculate energy and momentum 
if randDir  % randomly sample direction vectors 
    if isnumeric(p.Results.seed)
        v_hat = sampleDirVector(phi, theta, r_sphere, thresh, p.Results.seed); 
    else
        v_hat = sampleDirVector(phi, theta, r_sphere, thresh); 
    end
else
    v_hat = -r_0./vecnorm(r_0,2,2);  % initial momentum direction radial inward
end

m = 1.67262e-27;  % mass of proton [kg]
e = 1.6022e-19;  % charge on a proton, conversion from eV to J 
q = e;  % charge of particle [C]
c = 299792458;  % speed of light, m/s
B_0 = 1; 

R = m*c/q/B_0;  % Larmor radius
omega_0 = q*B_0/m; % cyclotron frequency
if ~isnumeric(KE)  % request random sampling kinetic energy
    KE = samplePowerlog(nRuns);  % default 1e6 eV to 1e9 eV
else
    KE = repmat(KE, [nRuns, 1]);  % convert to column vector
end
KE_J = KE*e; 
v = c*sqrt(1-(m*c^2./(KE_J+m*c^2)).^2);  % relativistic velocity to calculate tspan
p_hat = sqrt((1+KE_J./m/c^2).^2-1);
% p_hat = KE_J/m/c^2;  % KE as given by Paolo


if rel
    info = [r_0, v_hat.*p_hat];
else
    info = [r_0, v_hat.*v];  % NON-RELATIVISTIC EOM CHECK
end

%% Set global vars
GL('I', I);
GL('coil_mp', coil_mp);
GL('dL', dL);
if rel
    GL('r_0', R); 
else
    GL('r_0', 1); 
end

%% Begin ODE45 integrations
res = zeros(nRuns, 1); 

if plots
    streaks = cell(nRuns, 1);  % preallocate :)
end
for ii = 1:nRuns
    t = [0 2*r_sphere/v(ii)]; 
    s = omega_0*t; 
    
    IC = info(ii,:);
    if rel % ode45 integration
        IC(1:3) = IC(1:3)/R;  % scale position to non-dim. 
        [~, traj] = ode45(@eom_rad_rel, s, IC);
        traj = traj*R;  % re-scale position to dimensionalize
    else
        [~, traj] = ode45(@eom_rad, t, IC); 
    end
    trail = traj(:,1:3);  % xyz trail of points
    
    % check to see if the trail intersepts the sphere bounded by 'thresh'
    mags = vecnorm(trail,2,2); 
    [min_dist,ind] = min(mags);  % find index of nearest approach
    
    if ind == length(trail)  % nearest point to origin is at the end
        res(ii) = does_it_hit(trail(ind,:), trail(ind-1,:), thresh);
    elseif ind == 1
        res(ii) = does_it_hit(trail(ind,:), trail(ind+1,:), thresh); 
    else  % if nearest isn't at the end...
        if mags(ind+1)<mags(ind-1)  % check if next nearest point is ahead 
            res(ii) = does_it_hit(trail(ind,:), trail(ind+1,:), thresh);
        else  % nearest point must be behind
            res(ii) = does_it_hit(trail(ind,:), trail(ind-1,:), thresh);
        end
    end
    
    % store trajectories for plotting if true
    if plots
        streaks{ii} = trail; 
    end
end

defl_rate = (nRuns-sum(res))/nRuns; 

%% Finish plots
if plots
    f1 = figure();  
    LineSp = {'g:', 'r-'}; 
    LineTh = [0.5, 0.5];  % set up some linespec options in advance
    
    for ii = 1:nRuns  % plot trajectories
        x = streaks{ii}(:,1);
        y = streaks{ii}(:,2);
        z = streaks{ii}(:,3);
        plot3(x, y, z, LineSp{res(ii)+1}, 'LineWidth', LineTh(res(ii)+1)); 
    hold on;
    end
    
    r_plot = 25;
    plot_halbach(points, f1)
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