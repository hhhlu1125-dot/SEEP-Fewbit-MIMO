function [Sigma_q, Mu_q] = BEP_detection(RxSymbol_Q, Hreal, sigma2, sym, avg_symbol_energy, Q_bit, Delta, eta_b, iterNum)
        
        [nRow, nCol] = size(Hreal);
        TxAntNum = nCol / 2;
        RxAntNum = nRow / 2;

        delta=0.5; 

        % BEP检测器
        % BEP init
        Lambda = 0.0001*ones(2*TxAntNum, 1);     gamma = zeros(2*TxAntNum, 1);
     
        % Bussgang decomposition 
        [Hequal, Sigma_n] = bussgang_decomp(Hreal, sigma2, avg_symbol_energy, Q_bit, Delta, eta_b); 

        % Preprocessing
        Gram = Hequal'*(Sigma_n \ Hequal);
        Sigma_q = (Gram + diag(Lambda))\ eye(2*TxAntNum);  % 将MMSE结果作为预处理
        Mu_q = Sigma_q* (Hequal'*(Sigma_n \ RxSymbol_Q) + gamma);          % 

        for k = 1:iterNum 
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
    
            %EPD
            Sigma_q = (Gram + diag(Lambda))\ eye(2*TxAntNum);  % 将MMSE结果作为预处理
            Mu_q = Sigma_q* (Hequal'*(Sigma_n \ RxSymbol_Q) + gamma);          % 

        end