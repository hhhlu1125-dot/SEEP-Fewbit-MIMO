function [A, Sigma_n] = bussgang_decomp(r_channel, sigma2, avg_symbol_energy, b, Delta, eta_b)

    [nRow, nCol] = size(r_channel);
    MT2 = nCol;
    MR2 = nRow;

    % transmitted signal covariance
    Sigma_x = 0.5*avg_symbol_energy * eye(MT2);

    % covariance of unquantized received signal
    Sigma_r = r_channel * Sigma_x * r_channel.' + (sigma2 / 2) * eye(MR2);

    % numerical safeguard
    eps_val = 1e-12;
    d = diag(Sigma_r);
    d = max(d, eps_val);

    % Bussgang gain matrix V (diagonal)
    v_diag = zeros(MR2, 1);
    for k = 1:MR2
        tmp = 0;
        for i = 1:(2^b - 1)
            idx = i - 2^(b - 1);
            tmp = tmp + exp(-0.5 * Delta^2 * idx^2 / d(k));
        end
        v_diag(k) = (Delta / sqrt(2*pi)) * tmp / sqrt(d(k));
    end
    V = diag(v_diag);

    % equivalent channel
    A = V * r_channel;

    % equivalent noise covariance
    if b == 1
        D_inv = diag(1 ./ d);
        D_inv_sqrt = diag(1 ./ sqrt(d));

        term = D_inv_sqrt * Sigma_r * D_inv_sqrt;
        term = min(max(term, -1 + 1e-12), 1 - 1e-12);

        Sigma_n = (Delta^2 / (2*pi)) * ...
                  (asin(term) - term + (sigma2 / 2) * D_inv);
    else
        % approximate multi-bit covariance model
        Sigma_n = (sigma2 / 2) * (V * V.') + eta_b * diag(diag(Sigma_r));
        Sigma_n = Sigma_n + 1e-12 * eye(MR2);
    end
end