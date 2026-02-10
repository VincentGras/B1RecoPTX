function [b1p, b1m, residual] = xfl_proposed_reco(xfl_signal, alphaMat, betaMat, b1val)
% proposed voxel-wise reconstruction

% alphaMat: saturation encoding scheme, Volt.s
% betaMat: readout encoding scheme, Volt.s
% b1val: Initialization for the non-convex optimization of the sum-of-squares transmit profile

% output:
% residual: residual between the satTFL data and the signal calculated from the model using the estimated B1+ 
% b1p: Estimated B1+
% b1m: Estimated B1-

% b1p = x * u; u: phasor (unit-norm complex vector); x: sum-of-squares transmit profile (scalar)

    if nargin < 4
        b1val = [5, 10, 20]; % rad/s/volt 
    end

    sat = vecnorm(alphaMat, 2, 1) ~= 0; % cycles with saturation
    interf = ~sat; % cycles without saturation
    Nc = size(alphaMat, 1);
    assert(nnz(interf) >= Nc, 'expecting at least %d interf cycle', Nc);
    
    % estimate the phasor from the non-saturated cycles
    u = (betaMat(:, interf)') \ xfl_signal(interf);
    b1m = vecnorm(u, 2, 1);
    u = u / b1m;
    

    % estimate the sum-of-squares transmit profile from all satTFL data via nonlinear optimization
    alphaMat_sat = alphaMat(:, sat);
    oo = optimoptions('fmincon', 'Algorithm', 'interior-point', 'Display', 'off', 'SpecifyObjectiveGradient', true, 'Maxiter', 1000);
    
    linc_A = abs(alphaMat_sat' * u) / (pi); 
    linc_B = ones(nnz(sat), 1);

    fval = inf(size(b1val));
    x = zeros(size(b1val));

    for i = 1:numel(b1val)
        [x(i), fval(i)] = fmincon(@objective, b1val(i), linc_A, linc_B, [], [], 0.1, Inf, [], oo);
    end


    [~, i] = min(fval);
    x = x(i);

    
    b1p = x * u;
    b1m = b1m / x;
    
    Sfit = b1m * cos(abs(alphaMat' * b1p)) .* sin(betaMat' * b1p); % Signal computed from the model
    
    residual = norm(Sfit - xfl_signal)/ norm(xfl_signal);
    

    
   
    function [L, J] = objective(x)

        [S, JS] = signal_model(x);
        L = S - xfl_signal;
        X = abs(xfl_signal' * xfl_signal);
        L = abs(L' * L) / X;
        J = -2 * real(JS' * S) / X;
        
    end
        


    function [S, J] = signal_model(x)
        
        a = abs(alphaMat' * u);
        b = betaMat' * u;
        alpha = x * a;
        beta = x * b;
        
        sa = sin(alpha);
        ca = cos(alpha);

        sb = sin(beta);
        cb = cos(beta);

        
        
        s = sb .* ca;
        Js = b .* cb .* ca - a .* sb .* sa; 
        
        ss = s' * s;
        Jss = 2 * real(Js' * s); % + s' * Js;
        ssi = 1 / ss; 
        
        f = ((s' * xfl_signal) / (s' * s));
        Jf = ssi * (Js' * xfl_signal - ssi * Jss * (s' * xfl_signal));
        
        S = f  * s;
        J = Jf * s + f * Js;

        
    end


end
