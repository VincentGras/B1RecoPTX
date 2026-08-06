%% optimize nominal flip angles to minimize B1+ reconstruction error
clear all 
calcperf = @(val) mean(val(:, 1));  
calcperf_cell = @(cellval) cellfun(calcperf, cellval, 'UniformOutput', true);


db = '.\Avanti2_B1maps_Agar.map.mat'; % file with simulated B1+ fields (Nx x Ny x Nz x Ntx, T/V)
gyr = 42577000;
dataB1 = load(db);
data.M = dataB1.Mask;
data.b1 = map2vect(dataB1.Mask,  dataB1.B1map).' * (2*pi) * gyr; % (transform to vector form Ntx x Nvox, rad/s/V)

paramfile='.\Schemes\Sphere_vcc03_Acond2_Bcond2.xflparam.txt'; % file with encoding scheme 

rep = 10;  
noise_level = 1e-3;  %noise std


Vref = 165; % reference voltage (V)
A_ref = 40:10:150; % nominal saturation FA (deg)
B_ref = 2:2:30; % nominal readout FA  (deg)
[nom_alpha, nom_beta] = ndgrid(A_ref, B_ref);

result_reco = cell([size(nom_alpha)]);

for i = 1:numel(nom_alpha)

    result_reco{i} = perform_xfl_reco (paramfile, data, nom_alpha(i), nom_beta(i), Vref, @xfl_proposed_reco, rep, noise_level);
    fprintf('performance (alpha = %3d, beta = %3d): %5.2f \n', nom_alpha(i), nom_beta(i), calcperf(result_reco{i}));

end

% Plot the mean B1+ errors as a function of nominal FA pairs.
figure; 
imagesc(nom_alpha(:,1), nom_beta(1,:), calcperf_cell(squeeze(result_reco))');
caxis([0, 0.1])
axis xy

% Find the optimal alpha and beta values that minimize the mean B1+ reconstruction error
% subject to max(alpha) <= 180° and max(beta) <= 30°

perf_with_noise =  calcperf_cell(squeeze(result_reco));

[alphaMat, betaMat] = readXFLParametersFromFile(paramfile);

sat = alphaMat(:, any(alphaMat,1)) * Vref*1e-3/pi; % Volt.s/rad
read = betaMat * Vref*1e-3/pi;
FA_alpha_norm=abs(sat'*data.b1);
FA_beta_norm=abs(read'*data.b1);

a_valid = prctile(FA_alpha_norm(:), 99.9)*A_ref <= 180;
b_valid = prctile(FA_beta_norm(:), 99.9)*B_ref <=30;
valid_pairs = a_valid(:) & b_valid(:).';
perf = perf_with_noise(:);
perf(~valid_pairs) = Inf;
[~, idx] = min(perf(:));
alpha_optim = nom_alpha(idx); 
beta_optim = nom_beta(idx);


fprintf('%d %d\n', alpha_optim, beta_optim)


% write the optimized nominal FAs to the encoding scheme file
write_optim_angles = true;

if write_optim_angles

    lines = readlines(paramfile);

    newLines = [sprintf("#Vref = %.1f V", Vref); 
                sprintf("#alpha_nominal = %.1f deg", alpha_optim); 
                sprintf("#beta_nominal = %.1f deg", beta_optim)];
    
    start = ["#Vref ="; "#alpha_nominal ="; "#beta_nominal ="];
    

    missingLines = false(3, 1);
    
    for k = 1:3
        index = find(startsWith(strtrim(lines), start(k)), 1);
    
        if isempty(index)
            missingLines(k) = true;
        else
            lines(index) = newLines(k);
        end
    end
    
    index = find(startsWith(strtrim(lines), "#XFLparameters"), 1);

    if isempty(index)
        lines = ["#XFLparameters"; newLines(missingLines); lines];
    elseif any(missingLines)
        lines = [lines(1:index); newLines(missingLines); lines(index + 1:end)];
    end
    
    writelines(lines, paramfile);

end