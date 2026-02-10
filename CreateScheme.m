
addpath .\Functions\

folder_b1 = '.\Avanti2_B1maps_Agar.map.mat';  % file with simulated B1+ fields (Nx x Ny x Nz x Nch, T/V)
data = load(folder_b1);
B1map = data.B1map;
Mask = data.Mask;
gyr = 42577000;

B1_vec = map2vect(Mask, B1map).' * (2*pi) * gyr;  %  vector form Nch x Nvox, rad/s/V

par.Nsat=3; % number of presaturated acquisitions
par.maxcond_read = 2; % constraint on the condition number of the readout encoding matrix of non-saturated cycles (Y0)
par.maxcond_sat = 2; % constraint on the condition number of the saturation/readout encoding matrix of pre-saturated cycles (Y(sat)=X(sat))
par.FAsat_nom = 100; % nominal saturation flip angle (value in the protocol)
par.FAread_nom = 10; % nominal readout flip angle (value in the protocol)

run_vcc_optimize(B1_vec, par);