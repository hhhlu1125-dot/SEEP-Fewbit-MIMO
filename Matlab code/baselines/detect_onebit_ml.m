function x_hat = detect_onebit_ml(RxSymbol_Q, Hreal, sigma2, X_candidate_gpu)

    c = single(1.702);

    y = gpuArray(single(sign(RxSymbol_Q)));
    Hreal = gpuArray(single(Hreal));
    sigma2 = single(sigma2);

    X_real_gpu = [real(X_candidate_gpu);imag(X_candidate_gpu)];

    HX = Hreal * X_real_gpu;

    t = -c .* (y .* HX) ./ sqrt(sigma2 / 2);

    softplus = max(t, 0) + log1p(exp(-abs(t)));

    metric = sum(softplus, 1);

    [~, best_idx] = min(metric);

    x_hat = X_candidate_gpu(:, best_idx);

    x_hat = gather(x_hat);

end