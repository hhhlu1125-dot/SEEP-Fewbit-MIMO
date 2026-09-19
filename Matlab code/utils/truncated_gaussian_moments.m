function [mu_trunc, var_trunc] = truncated_gaussian_moments(mu, var, low, up)
%TRUNCATED_GAUSSIAN_MOMENTS
% Stable moments of truncated Gaussian for vector inputs.
%
% Input:
%   mu, var, low, up : vectors of the same size
%
% Output:
%   mu_trunc, var_trunc : mean and variance after truncation
    phi = @(x) exp(-0.5 .* x.^2) ./ sqrt(2*pi);
    Phi = @(x) 0.5 .* erfc(-x ./ sqrt(2));
    eps_var = 5e-7;
    eps_den = 1e-16;

    var = max(var, eps_var);
    sigma = sqrt(var);

    mu_trunc = zeros(size(mu));
    var_trunc = zeros(size(var));

    % ---------------------------------------------------------
    % 1) finite interval: (low, up]
    % ---------------------------------------------------------
    idx_finite = isfinite(low) & isfinite(up);

    if any(idx_finite)
        mu_f = mu(idx_finite);
        var_f = var(idx_finite);
        sigma_f = sigma(idx_finite);
        low_f = low(idx_finite);
        up_f  = up(idx_finite);

        a = (low_f - mu_f) ./ sigma_f;
        b = (up_f  - mu_f) ./ sigma_f;

        Phi_a = Phi(a);
        Phi_b = Phi(b);
        phi_a = phi(a);
        phi_b = phi(b);

        Z = Phi_b - Phi_a;
        Z = max(Z, eps_den);

        ratio1 = (phi_a - phi_b) ./ Z;
        ratio2 = (a .* phi_a - b .* phi_b) ./ Z;

        mu_tmp = mu_f + sigma_f .* ratio1;
        var_tmp = var_f .* (1 + ratio2 - ratio1.^2);

        mu_trunc(idx_finite) = mu_tmp;
        var_trunc(idx_finite) = max(var_tmp, eps_var);
    end

    % ---------------------------------------------------------
    % 2) lower truncation: (low, +inf)
    % ---------------------------------------------------------
    idx_lower = isfinite(low) & isinf(up) & (up > 0);

    if any(idx_lower)
        mu_l = mu(idx_lower);
        var_l = var(idx_lower);
        sigma_l = sigma(idx_lower);
        low_l = low(idx_lower);

        a = (low_l - mu_l) ./ sigma_l;

        denom = 1 - Phi(a);
        denom = max(denom, eps_den);

        lambda = phi(a) ./ denom;

        mu_tmp = mu_l + sigma_l .* lambda;
        var_tmp = var_l .* (1 + a .* lambda - lambda.^2);

        mu_trunc(idx_lower) = mu_tmp;
        var_trunc(idx_lower) = max(var_tmp, eps_var);
    end

    % ---------------------------------------------------------
    % 3) upper truncation: (-inf, up]
    % ---------------------------------------------------------
    idx_upper = isinf(low) & (low < 0) & isfinite(up);

    if any(idx_upper)
        mu_u = mu(idx_upper);
        var_u = var(idx_upper);
        sigma_u = sigma(idx_upper);
        up_u = up(idx_upper);

        b = (up_u - mu_u) ./ sigma_u;

        denom = Phi(b);
        denom = max(denom, eps_den);

        lambda = phi(b) ./ denom;

        mu_tmp = mu_u - sigma_u .* lambda;
        var_tmp = var_u .* (1 - b .* lambda - lambda.^2);

        mu_trunc(idx_upper) = mu_tmp;
        var_trunc(idx_upper) = max(var_tmp, eps_var);
    end

    % ---------------------------------------------------------
    % 4) fallback for abnormal cases
    % ---------------------------------------------------------
    idx_rest = ~(idx_finite | idx_lower | idx_upper);
    if any(idx_rest)
        mu_trunc(idx_rest) = mu(idx_rest);
        var_trunc(idx_rest) = max(var(idx_rest), eps_var);
    end

    % final safeguard
    bad_idx = ~isfinite(mu_trunc) | ~isfinite(var_trunc) | (var_trunc < eps_var);
    if any(bad_idx)
        mu_trunc(bad_idx) = min(max(mu(bad_idx), low(bad_idx)), up(bad_idx));
        var_trunc(bad_idx) = eps_var;
    end
end