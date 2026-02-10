addpath .\Functions\

paramfile='.\Schemes\InVivo_vcc03.xflparam.txt'; % file with encoding scheme
data = load('.\DataInVivo.mat').data; % acquired satTFL data (acquisition cycles x rx channels x voxels)
Mask = load('.\DataInVivo.mat').Mask; % Brain mask

% read the reference voltage and nominal flip angles alpha and beta from the text file to calculate RF pulse integrals IA and IB
fid = fopen(paramfile,'r');
numLines = 4;
for ii = 1:numLines
    line = fgetl(fid); 
    if (~isempty(regexpi(line, 'Vref')))
        a=regexp(line, '(\d*)', 'match', 'once' );
        Vref=str2num(a);
    end

    if (~isempty(regexpi(line, 'alpha_nominal')))
        a=regexp(line, '(\d*)', 'match', 'once' );
        alpha_nom=str2num(a);
    end

    if (~isempty(regexpi(line, 'beta_nominal')))
        a=regexp(line, '(\d*)', 'match', 'once' );
        beta_nom=str2num(a);
    end
end
fclose(fid);

% read saturation (alphaMat) and readout (betaMat) encoding matrices (unitless), Nch x Ncyc (Tx channels x acquisition cycles)
[alphaMat, betaMat] = readXFLParametersFromFile(paramfile);
 
RFPulseIntegral = beta_nom/180 * Vref * 1e-3; % Volt.s
RFPulseIntegralSat = alpha_nom/180 * Vref * 1e-3; % Volt.s
alphaMat = RFPulseIntegralSat * alphaMat; 
betaMat = RFPulseIntegral * betaMat;

Ntx = size(alphaMat, 1); % number of transmit channels
Nrx = size(data, 2); % number of receive channels
Nvox = size(data,3); % number of voxels
Ncycles = size(alphaMat, 2); % number of acquisition cycles
test_noise_thresh = 1 - Ncycles^-1;  % noise threshold


b1p = nan(Ntx, Nvox); % estimated B1+ maps
b1m = nan(Nrx, Nvox); % estimated B1- maps
residual = inf(1, Nvox); % residual between the satTFL data and the signal calculated from the model using the estimated B1+ 
b1val =  [5, 10, 20];% initialization for the non-convex optimization of the sum-of-squares transmit profile
t = nan(1, Nvox); % Test whether the voxel is dominated by noise

% reconstruction across all voxels
tic;
parfor i = 1:Nvox 
    [b1p(:,i), b1m(:,i), residual(i), t(i)] = xfl_Rx_reco(data(:,:,i), alphaMat, betaMat, b1val, test_noise_thresh);
    if (mod(i, 100) == 0 && residual(i) < Inf)
        fprintf('%d/%d, %f\n', i, Nvox, residual(i));
    end    
end
toc;

    
% normalize the phase
refphase = exp(1i * angle(sum(b1p, 1)));
b1p = b1p .* conj(refphase);
if (~isempty(b1m))
    b1m = b1m .* (refphase);
end

% visualization of the reconstructed B1+ maps

viewp(1, applyMask(Mask, vect2map(ones(size(Mask)), b1p)),'cscale',[0 30]); colormap('jet')

viewp(1, applyMask(Mask, vect2map(ones(size(Mask)), angle(b1p))),'cscale',[-pi pi]); colormap('jet')

