%% optimize nominal flip angles to minimize B1+ reconstruction error
clear all 
calcperf = @(val) mean(val(:, 1));  
calcperf_cell = @(cellval) cellfun(calcperf, cellval, 'UniformOutput', true);


db = '.\Avanti2_B1maps_Agar.map.mat'; % file with simulated B1+ fields (Nx x Ny x Nz x Nch, T/V)
gyr = 42577000;
dataB1 = load(db);
data.M = dataB1.Mask;
data.b1 = map2vect(dataB1.Mask,  dataB1.B1map).' * (2*pi) * gyr; % (transform to vector form Nch x Nvox, rad/s/V)

paramfile='.\Schemes\Sphere_vcc03_Acond2_Bcond2.xflparam.txt'; % file with encoding scheme

rep = 10;  % number of repetitions
noise_level = 1e-3;  %noise std


Vref=165; % reference voltage
A_ref = 40:10:150; % nominal saturation FA
B_ref = 2:2:30; % nominal readout FA
[nom_alpha, nom_beta] = ndgrid(A_ref, B_ref);

result_reco = cell([size(nom_alpha)]);

for i = 1:numel(nom_alpha)

    result_reco{i} = perform_xfl_reco (paramfile, ...
        data, nom_alpha(i), nom_beta(i), Vref, @xfl_proposed_reco, rep, noise_level);
    fprintf('simple reco : %20s(%5.1e, alpha = %3d, beta = %3d) : %5.2f \n', ...
            'avanti2_pairwise_modif_agar', noise_level, nom_alpha(i), nom_beta(i), ...
            calcperf(result_reco{i}));

end

% Plot the mean B1+ errors as a function of nominal FA pairs.
figure; 
imagesc(nom_alpha(:,1), nom_beta(1,:), calcperf_cell(squeeze(result_reco))');
caxis([0, 1])
axis xy

% Find the optimal alpha and beta that minimize the mean B1+ reconstruction error
perf_with_noise =  calcperf_cell(squeeze(result_reco));
[b1reco_perf, index_opt_alpha_beta] = min(perf_with_noise(:));
opt_alpha = nom_alpha(index_opt_alpha_beta);
opt_beta = nom_beta(index_opt_alpha_beta);
