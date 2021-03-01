function output = forFminSearch(x)

coilRadius = x(1);
halbachRadius = x(2);
AR = x(3);
geom = coil_racetrack(coilRadius, coilRadius/AR, 21);
[points, coil_mp, dL] = create_halbach(geom, 8, halbachRadius);
defl_rate = shielding_rate(points, coil_mp, dL);


% Fminsearch looks for minimum
output = 1-defl_rate;

end