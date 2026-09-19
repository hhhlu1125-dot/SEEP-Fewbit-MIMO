function xMean = detect_lmmse_qvb(RxSymbol_Q_c, RxSymbol_Q, HMat, const, prior, Q_bit, Delta, maxIter)
%DETECT_LMMSE_QVB
% LMMSE-QVB detector with known CSI.

[RxAntNum, TxAntNum] = size(HMat);
const = const(:).';
prior = prior(:).';
% Normalize prior
prior = prior / sum(prior);
eps_val = 1e-12;

% interval of RxSymbol
[r_low, r_up] = quantized_output_to_interval(RxSymbol_Q, Delta, Q_bit);
r_low_c = r_low(1:RxAntNum,1)+1j*r_low(RxAntNum+1:end,1);
r_up_c = r_up(1:RxAntNum,1)+1j*r_up(RxAntNum+1:end,1);

xMean = zeros(TxAntNum,1);

% tau_xi^1 = Var_p(x_i)[x_i]
priorMean = sum(prior .* const);
priorSecond = sum(prior .* abs(const).^2);
priorVar = real(priorSecond - abs(priorMean)^2);
priorVar = max(priorVar, eps_val);
xVar = repmat(priorVar, TxAntNum, 1);

% r_m^1 = y_m
rMean = RxSymbol_Q_c(:);
% tau_rm^1 = 0
rVar = zeros(RxAntNum,1);
% e = r - Hx
e = rMean - HMat * xMean;
eyeM = eye(RxAntNum);

% Store the z_i and precision used for the final MAP decision
zLast = zeros(TxAntNum,1);
precisionLast = zeros(TxAntNum,1);

%  VB iterations
for it = 1:maxIter

    Sigma_r = diag(max(real(rVar),0));
    Sigma_x = diag(max(real(xVar),0));
    Cpost = (norm(e)^2 / RxAntNum) * eyeM + Sigma_r + HMat * Sigma_x * HMat';
    % Remove tiny numerical non-Hermitian error
    Cpost = (Cpost + Cpost') / 2;
    Gamma = Cpost \ eyeM;
    Gamma = (Gamma + Gamma') / 2;

    for m = 1:RxAntNum
        oldRm = rMean(m);
        % [Gamma]_{m,m}
        precision_m = real(Gamma(m,m));
        precision_m = max(precision_m, eps_val);
        s_m = rMean(m)  - Gamma(m,:) * e / precision_m;
        [newRm, newRVar] = trunc_complex_moments(s_m, precision_m, r_low_c(m), r_up_c(m));
        rMean(m) = newRm;
        rVar(m) = newRVar;
        e(m) = e(m) - oldRm + newRm;

    end

    for i = 1:TxAntNum
        hi = HMat(:,i);
        % h_i^H Gamma h_i
        precision_i = real(hi' * Gamma * hi);
        precision_i = max(precision_i, eps_val);
        z_i = xMean(i) + hi' * Gamma * e / precision_i;
        [newXi, newXVar] = discrete_posterior_moments(z_i, precision_i, const, prior);
        oldXi = xMean(i);
        xMean(i) = newXi;
        xVar(i) = newXVar;
        e = e + hi * (oldXi - newXi);
        zLast(i) = z_i;
        precisionLast(i) = precision_i;

    end

end

end
