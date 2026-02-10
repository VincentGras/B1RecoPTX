function [b1p, b1m, residual, t] = xfl_Rx_reco(xfl_signal, alphaMat, betaMat, b1val, test_noise_thresh)
% proposed reconstruction for multi-Rx data

% xfl_signal: acquired satTFL data (acquisition cycles x rx channels x voxels)
% alphaMat: saturation encoding scheme, Volt.s
% betaMat: readout encoding scheme, Volt.s
% b1val: Initialization for the non-convex optimization of the sum-of-squares transmit profile
% test_noise_thresh: noise threshold

% output:
% residual: residual between the satTFL data and the signal calculated from the model using the estimated B1+ 
% b1p: Estimated B1+
% b1m: Estimated B1-
% t: test whether the voxel is dominated by noise

    Nrx = size(xfl_signal, 2); % number of receive channels
    Ntx = size(alphaMat, 1); % number of transmit channels
    b1p = NaN(Ntx, 1);
    b1m = NaN(Nrx, 1);
    residual = inf;

    if Nrx > 1

        [U, L, V] = svd(xfl_signal); % S = U*L*V': Singular value decomposition to determine the receive phasor (unit-norm vector)
        
        L = diag(L);
        t = test_noise(L);
        
        if t <= test_noise_thresh
            return;
        end
        
        % reco of B1+ and SoS receive profile 
        [b1p(:), b1m_, residual] = xfl_proposed_reco(U(:, 1), alphaMat, betaMat, b1val);  

        % B1- = SoS receive profile x receive phasor
        b1m(:) = (L(1, 1) * b1m_) * V(:, 1);
        
    else % Nrx == 1
        
        t = norm(xfl_signal);
        
        if t <= test_noise_thresh
            return;
        end
        

        [b1p(:), b1m(:), residual] = xfl_proposed_reco(xfl_signal, alphaMat, betaMat, b1val); 

     end

end


function t = test_noise(L)

    t = abs(L' * L);
    if t > 0
        t = abs(L(1))^2 / t;
    end
    
end