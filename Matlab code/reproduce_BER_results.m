%% BER VS SNR
clear;
clc;

rootDir = fileparts(mfilename('fullpath'));

addpath(fullfile(rootDir, 'detectors'));
addpath(fullfile(rootDir, 'baselines'));
addpath(fullfile(rootDir, 'utils'));

% profile on;
EbN0db      = 0:5:30;
SampleNum   = 1000000;       % max sample no.
isCorr = 0;
rho = 0.5;  % related channel parameter
max_frame   = 100;
ModType 	= 2;
codeRate = 1;
iterNum = 5;
iterNum_QVB = 10;
Q_bit = 1;
K = 0.3; % rician factor

% Detection Parameter Set
TxAntNum    = 8;
RxAntNum    = 32;
N 			= 1; 						% Block length
InfoLen 	= N * ModType * TxAntNum; 	% Information length

ReceiverType = 'ApproSEEP';

switch(ModType)
    case 2
        QAM=4;
        avg_symbol_energy = 2;
        sym=[-1 +1]';
    case 4
        QAM=16;
        avg_symbol_energy=10;
        sym = [-1 -3 +1 +3]';
    case 6
        QAM=64;
        avg_symbol_energy=42;
        sym=[-3 -1 -5 -7 +3 +1 +5 +7]';
    case 8
        QAM=256;
        avg_symbol_energy=170;
        sym=[-5.0 -7.0 -3.0 -1.0 -11.0 -9.0 -13.0 -15.0 ...
              5.0  7.0  3.0  1.0  11.0  9.0  13.0  15.0]';
end

[I, Q] = meshgrid(sym, sym);
sym_complex = I(:) + 1j * Q(:);
prior = ones(size(sym_complex)) / numel(sym_complex);
% ML need
% M = length(sym_complex);
% numCandidate = M^TxAntNum;
% idx = dec2base(0:numCandidate-1, M, TxAntNum) - '0' + 1;
% X_candidate = sym_complex(idx).';
% X_candidate_gpu = gpuArray(single(X_candidate));

slen = length(sym);
H           = zeros(2*RxAntNum,2*TxAntNum);  
RxSymbol    = zeros(2*RxAntNum,1); % Received

tic
for nEN =1:length(EbN0db)
    
    % initialization 
    SNR = EbN0db(nEN);
    SNR_indB = 10^(-SNR / 10)./ RxAntNum;

    error_frame = 0;
    error_bits = 0;
    
    TxBitsLLR = zeros(InfoLen, 1);
    
    loop = 0;
    while error_frame < max_frame && loop < SampleNum
        
        loop = loop + 1;

        TxBits = randi([0, 1], [InfoLen, 1]);
        TxSymbol=modulation(TxBits,ModType,sym);
        
        % Channel Matrix
        if isCorr
            HMat = rician_mimo_channel(RxAntNum, TxAntNum, K);
        else
            HMat = sqrt(0.5) *(randn(RxAntNum, TxAntNum) + sqrt(-1)*randn(RxAntNum, TxAntNum)); %Rayleigh channel
        end        
       
        % Received and sigma^2
        sigma2 = avg_symbol_energy*(norm(HMat, "fro")^2) * SNR_indB;
        RxSymbol_c=HMat*TxSymbol+sqrt(sigma2/2)*(randn(size(HMat*TxSymbol))+1i*randn(size(HMat*TxSymbol)));
        Nv=sigma2/avg_symbol_energy;


        if (mod(loop,1000) == 0)
            fprintf("\nNow Iter: %d\tNow SNR: %d\tNow Sigma: %f\tNow Error Frame: %d\tNow Error Bits: %d", ...
                loop, EbN0db(nEN), Nv, error_frame, error_bits);
        end

        RxSymbol(1:RxAntNum,1)=real(RxSymbol_c);
        RxSymbol(RxAntNum+1:end,1)=imag(RxSymbol_c);
      
        input_power = 0.5*(TxAntNum*avg_symbol_energy+sigma2);
        [RxSymbol_Q, Delta, eta_b] = few_bit_quantizer(RxSymbol, Q_bit, input_power);

        % Complex channel matrix to real
        Hreal = [real(HMat), -imag(HMat); imag(HMat), real(HMat)];

        switch ReceiverType
            case 'BMMSE'
                Mu_q = BMMSE_detection(RxSymbol_Q, Hreal, sigma2, avg_symbol_energy, Q_bit, Delta, eta_b);
                symLLR = Mu_q;
            case 'EP'
                [Sigma_q, Mu_q] = EP_detection(RxSymbol_Q, Hreal, sigma2/2, sym, iterNum);
                symLLR = Mu_q;
            case 'BEP'
                [Sigma_q, Mu_q] = BEP_detection(RxSymbol_Q, Hreal, sigma2/2, sym, Q_bit, Delta, eta_b, iterNum);
                symLLR = Mu_q;
            case 'SEEP'
                [Sigma_q, Mu_q] = SEEP_detection(RxSymbol_Q, Hreal, sigma2/2, sym, Q_bit, Delta, iterNum);
                symLLR = Mu_q;
            case 'ApproSEEP'
                [Sigma_q, Mu_q] = ApproSEEP_detection(RxSymbol_Q, Hreal, sigma2/2, sym, Q_bit, Delta, iterNum);
                symLLR = Mu_q;
            case 'MF_QVB'
                RxSymbol_Q_c = RxSymbol_Q(1:RxAntNum,1) + 1j*RxSymbol_Q(RxAntNum+1:end,1);
                xHard = detect_mf_qvb(RxSymbol_Q_c, RxSymbol_Q, HMat, sym_complex, prior, Q_bit, Delta, iterNum_QVB);
                symLLR = [real(xHard); imag(xHard)];
            case 'LMMSE_QVB'
                RxSymbol_Q_c = RxSymbol_Q(1:RxAntNum,1) + 1j*RxSymbol_Q(RxAntNum+1:end,1);
                xHard = detect_lmmse_qvb(RxSymbol_Q_c, RxSymbol_Q, HMat, sym_complex, prior, Q_bit, Delta, iterNum_QVB);
                symLLR = [real(xHard); imag(xHard)];
            case 'ML'
                xHard = detect_onebit_ml(RxSymbol_Q, Hreal, sigma2, X_candidate);
                symLLR = [real(xHard); imag(xHard)];
        end

        Tx_bitLLR = sym2bit_LLR(symLLR,ModType,TxAntNum,0);    
        Tx_Bits_est = Tx_bitLLR > 0;

        error_bits=error_bits+nnz(Tx_Bits_est-TxBits);
        
        if nnz(Tx_Bits_est-TxBits)
            
            error_frame = error_frame + 1;
            
        end

    end
    
    BER(nEN) = error_bits / TxAntNum / ModType / loop;
    FER(nEN) = error_frame / loop;
    
end

toc
figure(1);
semilogy(EbN0db,BER,'Marker','o','LineWidth',1.5);
xlim([min(EbN0db),max(EbN0db)])
ylim([1e-5,1]);
xlabel('SNR(dB)');
ylabel('BER');

