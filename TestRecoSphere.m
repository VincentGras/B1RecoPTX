%% check the method performance for the given encoding scheme using simulated B1+ maps (spherical phantom)
clear all
addpath .\Functions\

calcperf = @(val) mean(val(:, 1));  
calcperf_cell = @(cellval) cellfun(calcperf, cellval, 'UniformOutput', true);



db = '.\Avanti2_B1maps_Agar.map.mat'; % file with simulated B1+ fields (Nx x Ny x Nz x Ntx, T/V)
gyr = 42577000;
dataB1 = load(db);
data.M = dataB1.Mask; % mask
data.b1 = map2vect(dataB1.Mask,  dataB1.B1map).' * (2*pi) * gyr; % (transform to vector form Ntx x Nvox, rad/s/V)

paramfile='.\Schemes\Sphere_vcc03_Acond2_Bcond2.xflparam.txt'; % file with encoding scheme

rep = 1; % number of repetitions
noise_level = 1e-3; %noise std

% read the reference voltage and nominal flip angles alpha and beta from the text file to calculate RF pulse integrals IA and IB

fid = fopen(paramfile,'r');
numLines = 4;
for ii = 1:numLines

    line = fgetl(fid); 
    if ~ischar(line)
        break;
    end

    if contains(line,"Vref")
        a = regexp(line,'\d+\.?\d*','match','once');
        Vref = str2double(a);
    end

    if contains(line,"alpha_nominal")
        a = regexp(line,'\d+\.?\d*','match','once');
        alpha_nom = str2double(a);
    end

    if contains(line,"beta_nominal")
        a = regexp(line,'\d+\.?\d*','match','once');
        beta_nom = str2double(a);
    end

end
fclose(fid);

% perform a reco from a synthetic noisy data with the reconstruction method
[error, b1p_est, b1m_est]  = perform_xfl_reco(paramfile, data, alpha_nom, beta_nom, Vref, @xfl_proposed_reco, rep, noise_level); 

b1_error=error(:,1); % B1+ error with respect to the ground truth
fit_error=error(:,2); % residual between the synthetic satTFL data and the signal predicted by the model using the estimated B1+ 

% visualization of the estimated B1+ maps and SoS transmit profile
view3p(applyMask(data.M, vect2map(data.M, abs(b1p_est(:,1:size(b1p_est,2)/rep)))),'cscale',[0 50]); colormap('jet') % B1+ ampl
view3p(applyMask(data.M, vect2map(data.M, angle(b1p_est(:,1:size(b1p_est,2)/rep)))),'cscale',[-pi pi]); colormap('jet') % B1+ phase
view3p(applyMask(data.M, vect2map(data.M, vecnorm(b1p_est(:,1:size(b1p_est,2)/rep),2,1))),'cscale',[0 50]); colormap('jet') % SoS transmit profile

% distribution of the B1+ error with respect to the ground truth
figure, histogram(b1_error,[0:0.001:1],'Normalization','pdf');


