function xMean = detect_mf_qvb(RxSymbol_Q_c, RxSymbol_Q, HMat, const, prior, Q_bit, Delta, maxIter)
%DETECT_MF_QVB Matched-filter QVB detector, known CSI.

[RxAntNum, TxAntNum] = size(HMat); 
% interval of RxSymbol
[r_low, r_up] = quantized_output_to_interval(RxSymbol_Q, Delta, Q_bit);
r_low_c = r_low(1:RxAntNum,1)+1j*r_low(RxAntNum+1:end,1);
r_up_c = r_up(1:RxAntNum,1)+1j*r_up(RxAntNum+1:end,1);

xMean = zeros(TxAntNum, 1);
xVar = repmat(sum(prior .* abs(const).^2) - abs(sum(prior .* const)).^2, TxAntNum, 1);
rMean = RxSymbol_Q_c(:);
rVar = zeros(RxAntNum, 1);
g = sum(abs(HMat).^2, 1).';
e = rMean - HMat * xMean;

for it = 1:maxIter
    denom = norm(e)^2 + sum(rVar) + sum(g .* xVar);
    gamma = RxAntNum / max(real(denom), 1e-12);

    s = rMean - e;
    [rNew, rVar] = trunc_complex_moments(s, gamma, r_low_c(:), r_up_c(:));
    e = e - rMean + rNew;
    rMean = rNew;

    for i = 1:TxAntNum
        z = xMean(i) + HMat(:, i)' * e / max(g(i), 1e-12);
        [xm, xv] = discrete_posterior_moments(z, g(i) * gamma, const, prior);
        old = xMean(i);
        xMean(i) = xm;
        xVar(i) = xv;
        e = e + HMat(:, i) * (old - xMean(i));
    end
end
end
