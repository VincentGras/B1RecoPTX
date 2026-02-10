function [performance, b1p_est, b1m_est] = perform_xfl_reco(paramfile, b1file, alpha_nom, beta_nom, Vref, reco, rep, noiselevel)
% assessment of the B1+ reconstruction performance: B1+ post-reco from
% synthetic noisy data using the given method and encoding scheme

% paramfile: text file containing the encoding matrices
% b1file: file with "ground truth" B1+ maps
% reco: reconstruction algorithm
% rep: number of repetitions
% noiselevel: noise standard deviation

if nargin < 7
    rep = 1;
end

if nargin < 8
    noiselevel = 1e-3;
end



if isstr(b1file) 

    gyr = 42577000;
    dataB1 = load(b1file);
    data.M = dataB1.Mask;
    data.b1 = map2vect(dataB1.Mask,  dataB1.B1map).' * (2*pi) * gyr; % Nch x Nvox, rad/s/V

else
    data = b1file;
end


% read saturation (alphaMat) and readout (betaMat) encoding matrices (unitless), Nch x Ncyc (Tx channels x acquisition cycles)
[alphaMat, betaMat] = readXFLParametersFromFile(paramfile);

RFPulseIntegralAlpha = alpha_nom/180 * Vref * 1e-3; % Volt.s
alphaMat = RFPulseIntegralAlpha * alphaMat;
RFPulseIntegralBeta = beta_nom/180 * Vref * 1e-3; % Volt.s
betaMat = RFPulseIntegralBeta * betaMat;



b1 = data.b1;
if rep > 1
    b1 = repmat(b1, [1, rep]);
end

% synthesize noisy satTFL data and perform B1+ post-reconstruction
[fitting_residual, b1_estim_error, b1p_est, b1m_est] = simu_xfl_reco(alphaMat, betaMat, b1, noiselevel, reco);

% performance estimated as the B1+ error with respect to the ground truth (b1_estim_error) 
% and the fitting residual between the satTFL data and the signal calculated from the model using the estimated B1+ (fitting_residual)
performance = [b1_estim_error(:), fitting_residual(:)];

end
