function [xMean, xVar, xHard] = discrete_posterior_moments(z, gamma, const, prior)
%DISCRETE_POSTERIOR_MOMENTS Posterior mean/variance for z = x + CN(0,1/gamma).

eps_val = single(1e-12);
z = z(:);
gamma = real(gamma(:));

if numel(gamma)==1
    gamma = repmat(gamma,size(z));
end

const = const(:).';
prior = prior(:).';

dist2 = abs(z - const).^2;
logw = log(prior + eps_val) - gamma .* dist2;
logw = logw - max(logw,[],2);
w = exp(logw);
w = w ./ sum(w,2);
xMean = sum(w .* const,2);
secondMoment = sum(w .* abs(const).^2,2);
xVar = real(secondMoment - abs(xMean).^2);
xVar = max(xVar,eps_val);
[~,idx] = max(w,[],2);
xHard = const(idx).';
end
