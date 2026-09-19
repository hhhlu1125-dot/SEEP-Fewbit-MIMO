function [Sigma_q, Mu_q] = ApproSEEP_detection(RxSymbol_Q, Hreal, sigma2, sym, Q_bit, Delta, iterNum)

        [nRow, nCol] = size(Hreal);
        TxAntNum = nCol / 2;
        RxAntNum = nRow / 2;

        % interval of RxSymbol
        [r_low, r_up] = quantized_output_to_interval(RxSymbol_Q, Delta, Q_bit);

        delta = 0.5; % damping factor
        % EP init
        % x-side Gaussian sites
        Lambda = 0.0001*ones(2*TxAntNum, 1);     gamma = zeros(2*TxAntNum, 1);

        % r-side Gaussian sites
        Beta = 0.0001*ones(2*RxAntNum, 1);     alpha = zeros(2*RxAntNum, 1);
    
        % Preprocessing
        M = 1 / sigma2 + Beta;
        M_inv = diag(1 ./ M);
        A_r = M_inv * Hreal / sigma2;
        b_r = M_inv * alpha;
        Sigma_q = (Hreal' * diag(Beta) * A_r + diag(Lambda)) \ eye(2*TxAntNum);  % 将MMSE结果作为预处理
        Mu_q = Sigma_q * (Hreal' * b_r / sigma2 + gamma);          % 
        Sigma_q_r = M_inv + diag(sum((A_r * Sigma_q) .* A_r, 2)); %仅计算对角线元素
        Mu_q_r = A_r * Mu_q + b_r;

        for k = 1:iterNum
            % refine posterio probability of unquantized RX symbol
            diagSigma_q_r = diag(Sigma_q_r);
            m2 = diagSigma_q_r ./ (1 - diagSigma_q_r .* Beta);
            theta = m2 .* (Mu_q_r ./ diagSigma_q_r - alpha);

            [mu_p_r, sigma2_p_r] = truncated_gaussian_moments(theta, m2, r_low, r_up);

            tempBeta = 1 ./ sigma2_p_r - 1 ./ m2;
            tempalpha = mu_p_r ./ sigma2_p_r - theta ./ m2;
            Beta_new = Beta;
            alpha_new = alpha;      
            valid_idx = tempBeta > 5e-7;
            Beta_new(valid_idx) = tempBeta(valid_idx);
            alpha_new(valid_idx) = tempalpha(valid_idx);
            Beta = delta* Beta_new + (1-delta)* Beta;
            alpha = delta* alpha_new + (1-delta)* alpha;
            
            % refine posterio probability of TX symbol
            diagSigma = diag(Sigma_q);
            h2 = diagSigma ./ (1 - diagSigma .* Lambda);
            t  = h2 .* (Mu_q ./ diagSigma - gamma);
            
            [mu_p, sigma2_p] = prior_moments(t, h2, sym);

            % site update
            tempLamda = 1 ./ sigma2_p - 1 ./ h2;
            tempgamma = mu_p ./ sigma2_p - t ./ h2;
            Lambda_new = Lambda;
            gamma_new = gamma;            
            valid_idx = tempLamda > 5e-7;
            Lambda_new(valid_idx) = tempLamda(valid_idx);
            gamma_new(valid_idx) = tempgamma(valid_idx);
            Lambda = delta* Lambda_new + (1-delta)* Lambda;
            gamma = delta* gamma_new + (1-delta)* gamma;

            M = 1 / sigma2 + Beta;
            M_inv = diag(1 ./ M);
            Gram = Hreal'*diag(Beta)*M_inv*Hreal/sigma2;
            ytilde = Hreal'*M_inv*alpha/sigma2;
            Sigma_q = (Gram + diag(Lambda))\ eye(2*TxAntNum);  
            Mu_q = Sigma_q * (ytilde + gamma);         
    
            A_r = M_inv * Hreal / sigma2;
            b_r = M_inv * alpha;
            Sigma_q_r = M_inv + diag(sum((A_r * Sigma_q) .* A_r, 2)); %仅计算对角线元素
            Mu_q_r = A_r*Mu_q + b_r;
        end