%% Current vs Protection Script
clear

%% FMINSEARCH
thresh = .4;
for ii = 1:length(thresh)
    x{ii} = fminsearchbnd(@current_protection,[1, 1.5, 3, 1e6, thresh(ii)], [0.6, 0.5, 2, 1e6-1, thresh(ii)-1e-20],[5, 3, 10,  1e6+1, thresh(ii)+1e-20]);
end

AR = x{1}(1);
major_r = x{1}(2);
Hradius = x{1}(3);
I = x{1}(4);
thresh = x{1}(5);
geom = coil_racetrack(major_r, major_r/AR, 21);
[points, coil_mp, dL] = create_halbach(geom, 8, Hradius);
plot_halbach(points)
shielding = rel_shielding_rate(points, coil_mp, dL, I)
