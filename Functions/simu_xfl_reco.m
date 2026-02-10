function [fitting_residual, b1_estim_error, b1p_est, b1m_est] = simu_xfl_reco(alphaMat, betaMat, b1, noise_level, reco)
% assessment of the B1+ reconstruction performance: B1+ post-reco from synthetic noisy data using the given method

% alphaMat: saturation encoding scheme, Volt.s
% betaMat: readout encoding scheme, Volt.s
% b1: "ground truth" B1+ maps, rad/s/V
% noise_level: noise standard deviation
% reco: reconstruction algorithm

% output:
% fitting_residual: residual between the satTFL data and the signal calculated from the model using the estimated B1+ 
% b1_estim_error: the B1+ error with respect to the ground truth
% b1p_est: Estimated B1+ maps
% b1m_est: Estimated B1- profile

if (nargin < 5)
    reco = @xfl_proposed_reco;
end


alpha = alphaMat' * b1; % saturation FA, Ncyc x Nvox, rad
beta = betaMat' * b1; % readout FA, Ncyc x Nvox, rad

S0 = cos(abs(alpha)) .* sin(beta); % signal model, Ncyc x Nvox

cmplxnoise = @(noiselevel, sz) (randn(sz) + 1i * randn(sz)) * (noiselevel/sqrt(2));

S = S0 + cmplxnoise(noise_level, size(S0)); % noisy signal, Ncyc x Nvox

fitting_residual = inf(1, size(S, 2));
b1_estim_error = inf(1, size(S, 2));
b1p_est = nan(size(alphaMat, 1), size(S, 2));
b1m_est= nan(1, size(S, 2));

% reconstruction across all voxels
parfor i = 1:size(S, 2) 

    [b1p_est_, b1m_est_, fitting_residual(i)] = reco(S(:,i), alphaMat, betaMat);
    ph = angle(b1(:,i)'*b1p_est_);
    b1_estim_error(i) = norm(exp(1i * ph) * b1(:,i) -  b1p_est_) / norm(b1(:,i));
    b1p_est(:, i) = b1p_est_;
    if mod(i, 1000) == 0
        fprintf('%d / %d, fitting_residual = %f\n', i, size(S, 2), fitting_residual(i));
    end

    b1m_est(i) = b1m_est_;
end







