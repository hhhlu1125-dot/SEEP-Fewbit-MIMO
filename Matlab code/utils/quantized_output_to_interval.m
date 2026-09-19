function [r_low, r_up] = quantized_output_to_interval(RxSymbol_Q, Delta, Q_bit)
        tau_low  = (-2^(Q_bit-1) + 1) .* Delta;
        tau_high = ( 2^(Q_bit-1) - 1) .* Delta;
        r_up  = inf(size(RxSymbol_Q));
        r_low = -inf(size(RxSymbol_Q));
        idx_up  = RxSymbol_Q < tau_high;
        idx_low = RxSymbol_Q > tau_low;
        r_up(idx_up)   = RxSymbol_Q(idx_up)   + Delta/2;
        r_low(idx_low) = RxSymbol_Q(idx_low)  - Delta/2;
end