function [mu_p, sigma2_p] = prior_moments(mu, var, sym)

            % avoid numerical issues
            var = max(var, 1e-12);
            
            % log-probabilities for all dimensions and all constellation points
            logprob = - (sym.' - mu).^2 ./ (2 * var);
            logprob = logprob - max(logprob, [], 2);
            
            prob = exp(logprob);
            prob = prob ./ sum(prob, 2);
            
            % tilted moments
            mu_p = prob * sym;
            sigma2_p = prob * (sym.^2) - mu_p.^2;
            sigma2_p = max(sigma2_p, 5e-7);
end
