function [b1p, b1m, residual] = xfl_simple_reco(xfl_signal, alphaMat, betaMat, b1val)

sat = vecnorm(alphaMat, 2, 1) ~= 0;
assert(any(sat), 'expecting at least one saturated cycle');
assert(nnz(sat) == 1, 'expecting only one saturated cycle');
assert(sat(end), 'saturated image must be the last one');


alpha_ref = acos(real(xfl_signal(end) / xfl_signal(end - 1))) ;
beta_ref = alpha_ref * (norm(betaMat(:, end))/norm(alphaMat(:, end)));
beta_interf = xfl_signal(1:end-2)/xfl_signal(end-1) * beta_ref;
b1p = (betaMat(:, 1:end-2)') \ beta_interf;

Sfit = sin(betaMat' * b1p) .* cos(abs(alphaMat' * b1p));
b1m = (Sfit' * xfl_signal) / (Sfit' * Sfit); 
Sfit = b1m * Sfit;
residual = norm(Sfit - xfl_signal) / norm(xfl_signal);


