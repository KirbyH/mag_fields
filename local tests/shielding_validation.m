% simple code to check if/when shielding values approximately converge

r_maj = 2; 
AR = 1; 
r_H = 5; 
geom = coil_racetrack(2, 1, 51); 
[points, coil_mp, dL] = create_halbach(geom, 8, 5); 
I = 5e6; 
KE = 1e8;  % using powerlog

N = 3:12; 
N = 2.^N; 
nTests = length(N); 
eff = zeros(nTests, 1); 

for ii = 1:nTests
    eff(ii) = rand_shielding(points, coil_mp, dL, ...
        'I', I, 'N', N(ii), 'thresh', 1.5, 'KE', 'powerlog'); 
end

figure; 
semilogx(N, (1-eff)*100); 
grid on; 
xlabel('Number of particles thrown [-]'); 
ylabel('Transmission probability [\%]'); 
ylim([0 100]); 
legend(sprintf('$r_{maj}=$ %.1f, $r_H=$ %.1f, $AR=$ %.1f ...\n$KE=$ %s, $I=$ %.1s', ...
    r_maj, r_H, AR, KE, I));
