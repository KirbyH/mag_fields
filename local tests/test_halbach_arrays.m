% clc; clear; close all; 

figure; 
geom = coil_geom(1, 0.5, 40); 
[points, coil_mp, dL] = create_halbach(geom, 1, 0); 
plot_halbach(points); 