%% check the method performance for the given encoding scheme using simulated B1+ maps (spherical phantom)

addpath .\Functions\

calcperf = @(val) mean(val(:, 1));  
calcperf_cell = @(cellval) cellfun(calcperf, cellval, 'UniformOutput', true);



db = '.\Avanti2_B1maps_Agar.map.mat'; % file with simulated B1+ fields (Nx x Ny x Nz x Nch, T/V)
gyr = 42577000;
dataB1 = load(db);
data.M = dataB1.Mask; % mask
data.b1 = map2vect(dataB1.Mask,  dataB1.B1map).' * (2*pi) * gyr; % (transform to vector form Nch x Nvox, rad/s/V)

paramfile='.\Schemes\Sphere_vcc03_Acond2_Bcond2.xflparam.txt'; % file with encoding scheme

rep = 1; % number of repetitions
noise_level = 1e-3; %noise std

% read the reference voltage and nominal flip angles alpha and beta from the text file to calculate RF pulse integrals IA and IB

fid = fopen(paramfile,'r');
numLines = 4;
for ii = 1:numLines
    line = fgetl(fid); 
    if (~isempty(regexpi(line, 'Vref')))
        a=regexp(line, '(\d*)', 'match', 'once' );
        Vref=str2num(a); % V
    end

    if (~isempty(regexpi(line, 'alpha_nominal')))
        a=regexp(line, '(\d*)', 'match', 'once' );
        alpha_nom=str2num(a); % degree
    end

    if (~isempty(regexpi(line, 'beta_nominal')))
        a=regexp(line, '(\d*)', 'match', 'once' );
        beta_nom=str2num(a); % degree
    end
end
fclose(fid);

% perform a reco from a synthetic noisy data with the reconstruction method
[error, b1p_est, b1m_est]  = perform_xfl_reco(paramfile, data, alpha_nom, beta_nom, Vref, @xfl_proposed_reco, rep, noise_level);

b1_error=error(:,1); % B1+ error with respect to the ground truth
fit_error=error(:,2); % residual between the satTFL data and the signal calculated from the model using the estimated B1+ 

% visualization of the estimated B1+ maps and SoS transmit profile
view3p(applyMask(data.M, vect2map(data.M, b1p_est(:,1:size(b1p_est,2)/rep))),'cscale',[0 50]); colormap('jet')
view3p(applyMask(data.M, vect2map(data.M, vecnorm(b1p_est(:,1:size(b1p_est,2)/rep),2,1))),'cscale',[0 50]); colormap('jet')

% distribution of the B1+ error with respect to the ground truth
figure, histogram(b1_error,[0:0.001:1],'Normalization','pdf');


