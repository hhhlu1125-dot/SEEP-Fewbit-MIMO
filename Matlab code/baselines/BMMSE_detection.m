function Mu_q = BMMSE_detection(RxSymbol_Q, Hreal, sigma2, avg_symbol_energy, Q_bit, Delta, eta_b)
        
        [nRow, nCol] = size(Hreal);
        TxAntNum = nCol / 2;
        RxAntNum = nRow / 2;

        % Bussgang decomposition 
        [Hequal, Sigma_n] = bussgang_decomp(Hreal, sigma2, avg_symbol_energy, Q_bit, Delta, eta_b); 
        % BMMSE 
        Sigma_x = (avg_symbol_energy/2) * eye(2*TxAntNum);
        B = Hequal * Sigma_x * Hequal' + Sigma_n;
        Mu_q = Sigma_x * Hequal' * (B \ RxSymbol_Q);
        end