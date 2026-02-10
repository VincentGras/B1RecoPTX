function [b1p, b1m, residual] = xfl_pairwisecp_reco(xfl_signal, alphaMat, betaMat)

Nc = size(alphaMat, 1);

sat = vecnorm(alphaMat, 2, 1) ~= 0;
interf = ~sat;
assert(size(alphaMat, 2) == 2 * Nc + 1, 'expecting %d cycles', 2*Nc + 1)
assert(nnz(sat) == Nc, 'expecting %d one saturated cycle', Nc);

% assume that the last interf is cp
cp = find(interf, 1, 'last');
interf(cp) = false;


r = real(xfl_signal(sat) ./ xfl_signal(cp));
r(r > 1) = 1;
r(r < -1) = -1;
alpha_interf = acos(r) .* exp(1i * angle(xfl_signal(interf))); 
beta_interf = alpha_interf .* ...
    (vecnorm(betaMat(:, interf), 2, 1) ./ vecnorm(alphaMat(:, sat), 2, 1))';
b1p = betaMat(:, interf)' \ beta_interf;

Sfit = sin(betaMat' * b1p) .* cos(abs(alphaMat' * b1p));
b1m = (Sfit' * xfl_signal) / (Sfit' * Sfit); 
Sfit = b1m * Sfit;
residual = norm(Sfit - xfl_signal) / norm(xfl_signal);
        