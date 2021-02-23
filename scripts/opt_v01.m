%% Plot Error Space for Optimizations
% Matt Tuman
% 2/22/2021
% This script utilizes the numerous functions graciously constructed by Mr.
% Kirby Heck to plot the objective function space for varying variables
clear all
close all

% How fine the grid is
density = 2^6;

%% Investigate Coil Geometry 
%%{
% Pre-allocates memmory :)
Objective = zeros(density,density);
% Limits for the horizontal and vertical radius
r = linspace(0.1, 20, density);
% Determines the number of points on the coil geometry
n_points = 2^4+1;
% Deterimens the Halbach radius
halbach_r = 10;
% Deterimes the number of coils
n_coils = 8;
% Determines the current
I = 1e6;

for ii = 1:density
    for jj = 1:density
        clear ll
        % Calculate the coil geometry
        geom = coil_geom(r(ii),r(jj),n_points);
        % Calculate the Halbach Array Geometry
        [points, coil_mp, dL] = create_halbach(geom, n_coils, halbach_r); 
        % calculate the coil forces on each point in the Halbach Array
        forces = coil_forces(coil_mp, dL, I); 
        % Loop through each coil
        for i = 1:n_coils
            % Index to grab the correct coil points
            span = 1+(n_points-1)*(i-1) : (n_points-1)*(i);
            % Sum up the forces for each coil
            sorry_kirby = sum(forces(span,:));
            % Compute the euclidean norm of the force on each coil
            ll(i) = norm(sorry_kirby);
        end
        % Sum up all of the euclidean norms on the coils for a specific
        % coil geometry and set it to a log scale
        Objective(jj,ii) = log(sum(ll));
    end
end

% Plot the results
figure
pcolor(r,r,Objective);
title('Force Norm Space Varying Coil Geometry')
xlabel('Horizontal Radius [m]')
ylabel('Vertical Radius [m]')
colorbar
%}

%% Investigate Horizontal Coil Radius and Halbach Geometry
%%{
clearvars -except density
% Pre-allocates memmory :)
Objective = zeros(density,density);
% Limits for the horizontal and vertical radius
r = linspace(0.1, 10, density);
% Kept
r_set = 1;
% Determines the number of points on the coil geometry
n_points = 2^4+1;
% Deterimens the Halbach radius
halbach_r = linspace(2,10, density);
% Deterimes the number of coils
n_coils = 8;
% Determines the current
I = 1e6;

for ii = 1:density
    for jj = 1:density
        clear ll
        % Calculate the coil geometry
        geom = coil_geom(r(ii),r_set,n_points);
        % Calculate the Halbach Array Geometry
        [points, coil_mp, dL] = create_halbach(geom, n_coils, halbach_r(jj)); 
        % calculate the coil forces on each point in the Halbach Array
        forces = coil_forces(coil_mp, dL, I); 
        % Loop through each coil
        for i = 1:n_coils
            % Index to grab the correct coil points
            span = 1+(n_points-1)*(i-1) : (n_points-1)*(i);
            % Sum up the forces for each coil
            sorry_kirby = sum(forces(span,:));
            % Compute the euclidean norm of the force on each coil
            ll(i) = norm(sorry_kirby);
        end
        % Sum up all of the euclidean norms on the coils for a specific
        % coil geometry and set it to a log scale
        Objective(jj,ii) = log(sum(ll));
    end
end

% Plot the results
figure
pcolor(r,halbach_r,Objective);
title('Force Norm Space Varying Horizontal Coil Radius and Halbach Array Radius')
xlabel('Horizontal Radius [m]')
ylabel('Halbach Radius [m]')
colorbar
%}

%% Investigate Vertical Coil Radius and Halbach Geometry

clearvars -except density
% Pre-allocates memmory :)
Objective = zeros(density,density);
% Limits for the horizontal and vertical radius
r = linspace(0.1, 10, density);
% Kept
r_set = 1;
% Determines the number of points on the coil geometry
n_points = 2^4+1;
% Deterimens the Halbach radius
halbach_r = linspace(1,20, density);
% Deterimes the number of coils
n_coils = 8;
% Determines the current
I = 1e6;

for ii = 1:density
    for jj = 1:density
        clear ll
        % Calculate the coil geometry
        geom = coil_geom(r_set,r(ii),n_points);
        % Calculate the Halbach Array Geometry
        [points, coil_mp, dL] = create_halbach(geom, n_coils, halbach_r(jj)); 
        % calculate the coil forces on each point in the Halbach Array
        forces = coil_forces(coil_mp, dL, I); 
        % Loop through each coil
        for i = 1:n_coils
            % Index to grab the correct coil points
            span = 1+(n_points-1)*(i-1) : (n_points-1)*(i);
            % Sum up the forces for each coil
            sorry_kirby = sum(forces(span,:));
            % Compute the euclidean norm of the force on each coil
            ll(i) = norm(sorry_kirby);
        end
        % Sum up all of the euclidean norms on the coils for a specific
        % coil geometry and set it to a log scale
        Objective(jj,ii) = log(sum(ll));
    end
end

% Plot the results
figure
pcolor(r,halbach_r,Objective);
title('Force Norm Space Varying Vertical Coil Radius and Halbach Array Radius')
xlabel('Vertical Radius [m]')
ylabel('Halbach Radius [m]')
colorbar