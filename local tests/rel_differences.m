m = 1.67262e-27;  % mass of proton [kg]
e = 1.6022e-19;  % charge on a proton, conversion from eV to J 
c = 299792458;  % speed of light, m/s

KE_eV = (linspace(3, 11, 101));  % energy, eV
KE_eV = 10.^KE_eV; 
KE_J = KE_eV*e; 

%% non-relativistic vs relativistic velocity plots
v_nr = sqrt(2*KE_J/m); 
v_rel = c*sqrt(1-(m*c^2./(KE_J+m*c^2)).^2); 

figure; 
loglog(KE_eV, v_nr, 'b--', KE_eV, v_rel, 'rx'); hold on; 
loglog(xlim, [c c], 'k:'); 
grid on; 
xlabel('Energy [eV]'); 
ylabel('Speed [m/s]'); 

%% non-relativistic vs relativistic momentum
p_nr = m*v_nr; 
p_hat = sqrt((1+KE_J/m/c^2).^2-1);
p_rel = p_hat*m*c; 

figure; 
loglog(KE_eV, p_nr, 'b--', KE_eV, p_rel, 'rx'); 
grid on; 
xlabel('Energy [eV]'); 
ylabel('Momentum [N-s]'); 