function [m, v] = trunc_complex_moments(mu, gamma, low, up)
%TRUNC_COMPLEX_MOMENTS Moments of CN(mu, 1/gamma) truncated by rectangle.
% The real and imaginary components are independent N(mu_component, 1/(2gamma)).
eps_val = single(1e-12);
gamma = max(real(gamma), eps_val);
[mr, vr] = trunc_real_moments(real(mu), gamma, real(low), real(up));
[mi, vi] = trunc_real_moments(imag(mu), gamma, imag(low), imag(up));

m = mr + 1j * mi;
v = max(real(vr+vi),eps_val);
end

function [m, v] = trunc_real_moments(mu, gamma, low, up)
eps_val = single(1e-12);
sigma = sqrt(1 ./ (2 * gamma));
alpha = (low - mu) ./ sigma;
beta = (up - mu) ./ sigma;

PhiA = normcdf_stable(alpha);
PhiB = normcdf_stable(beta);
Z = max(PhiB - PhiA, eps_val);

phiA = normpdf_stable(alpha);
phiB = normpdf_stable(beta);

lambda = (phiA - phiB) ./ Z;
m = mu + sigma .* lambda;

alphaPhi = alpha .* phiA;
betaPhi = beta .* phiB;
alphaPhi(~isfinite(alphaPhi)) = 0;
betaPhi(~isfinite(betaPhi)) = 0;

v = sigma.^2 .* (1 + (alphaPhi - betaPhi) ./ Z - lambda.^2);
v = max(real(v), eps_val);
end

function y = normpdf_stable(x)
x = real(x);
y = exp(-0.5.*x.^2) ./ sqrt(single(2*pi));
y(~isfinite(x))=0;
end

function y = normcdf_stable(x)
x = real(x);
y = 0.5 .* (1 + erf(x./sqrt(single(2))));
y(x==inf)=1;
y(x==-inf)=0;
end
