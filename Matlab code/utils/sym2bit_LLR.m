function TxBitsLLR = sym2bit_LLR(symLLR,ModType,TxAntNum, isSoft)

TxBitsLLR = zeros(TxAntNum*ModType,1);

if isSoft

    switch ModType
    
        case 4 % 16-QAM
    
            for i = 1:TxAntNum
                TxBitsLLR((i-1)*ModType+1) = max(symLLR(3:4,i)) - max(symLLR(1:2,i)); % 1/0 LLR；
                TxBitsLLR((i-1)*ModType+2) = max(symLLR([2 4],i)) - max(symLLR([1 3],i)); % 1/0 LLR；
                TxBitsLLR((i-1)*ModType+3) = max(symLLR(3:4,i+TxAntNum)) - max(symLLR(1:2,i+TxAntNum)); % 1/0 LLR；
                TxBitsLLR((i-1)*ModType+4) = max(symLLR([2 4],i+TxAntNum)) - max(symLLR([1 3],i+TxAntNum)); % 1/0 LLR；
            end 
    
        case 6 % 64-QAM
    
            for i = 1:TxAntNum
    
    
                TxBitsLLR((i-1)*ModType+1) = max(symLLR(5:8,i)) - max(symLLR(1:4,i));
                TxBitsLLR((i-1)*ModType+2) = max(symLLR([3 4 7 8],i)) - max(symLLR([1 2 5 6],i));
                TxBitsLLR((i-1)*ModType+3) = max(symLLR([2 4 6 8],i)) - max(symLLR([1 3 5 7],i));
    
                TxBitsLLR((i-1)*ModType+4) = max(symLLR(5:8,i+TxAntNum)) - max(symLLR(1:4,i+TxAntNum));
                TxBitsLLR((i-1)*ModType+5) = max(symLLR([3 4 7 8],i+TxAntNum)) - max(symLLR([1 2 5 6],i+TxAntNum));
                TxBitsLLR((i-1)*ModType+6) = max(symLLR([2 4 6 8],i+TxAntNum)) - max(symLLR([1 3 5 7],i+TxAntNum));
    
            end
    
        case 8 % 256QAM
    
            for i = 1:TxAntNum
    
                TxBitsLLR((i-1)*ModType+1) = max(symLLR(9:16,i)) - max(symLLR(1:8,i));
                TxBitsLLR((i-1)*ModType+2) = max(symLLR([5:8 13:16],i)) - max(symLLR([1:4 9:12],i));
                TxBitsLLR((i-1)*ModType+3) = max(symLLR([3 4 7 8 11 12 15 16],i)) - max(symLLR([1 2 5 6 9 10 13 14],i));
                TxBitsLLR((i-1)*ModType+4) = max(symLLR(2:2:16,i)) - max(symLLR(1:2:16,i));
    
                TxBitsLLR((i-1)*ModType+5) = max(symLLR(9:16,i+TxAntNum)) - max(symLLR(1:8,i+TxAntNum));
                TxBitsLLR((i-1)*ModType+6) = max(symLLR([5:8 13:16],i+TxAntNum)) - max(symLLR([1:4 9:12],i+TxAntNum));
                TxBitsLLR((i-1)*ModType+7) = max(symLLR([3 4 7 8 11 12 15 16],i+TxAntNum)) - max(symLLR([1 2 5 6 9 10 13 14],i+TxAntNum));
                TxBitsLLR((i-1)*ModType+8) = max(symLLR(2:2:16,i+TxAntNum)) - max(symLLR(1:2:16,i+TxAntNum));
            end

    end

else


    switch ModType

        case 2 % QPSK
            for i = 1:TxAntNum
                TxBitsLLR((i-1)*ModType+1) = symLLR(i);
                TxBitsLLR((i-1)*ModType+2) = symLLR(i+TxAntNum);
            end
    
        case 4 % 16-QAM
    
            for i = 1:TxAntNum
                TxBitsLLR((i-1)*ModType+1) = symLLR(i); % 1/0 LLR；
                TxBitsLLR((i-1)*ModType+2) = abs(real(symLLR(i))) - 2; % 1/0 LLR；2
                TxBitsLLR((i-1)*ModType+3) = symLLR(i+TxAntNum); % 1/0 LLR；
                TxBitsLLR((i-1)*ModType+4) = abs(symLLR(i+TxAntNum)) - 2; % 1/0 LLR；2
            end
    
        case 6 % 64-QAM
    
            for i = 1:TxAntNum
    
    
                TxBitsLLR((i-1)*ModType+1) = symLLR(i);
                TxBitsLLR((i-1)*ModType+2) = abs(TxBitsLLR((i-1)*ModType+1)) - 4; %4
                TxBitsLLR((i-1)*ModType+3) = abs(TxBitsLLR((i-1)*ModType+2)) - 2; %2
    
                TxBitsLLR((i-1)*ModType+4) = symLLR(i+TxAntNum);
                TxBitsLLR((i-1)*ModType+5) = abs(TxBitsLLR((i-1)*ModType+4)) - 4; %4
                TxBitsLLR((i-1)*ModType+6) = abs(TxBitsLLR((i-1)*ModType+5)) - 2; %2
    
            end
    
        case 8 % 256QAM
    
            for i = 1:TxAntNum
    
                TxBitsLLR((i-1)*ModType+1) = symLLR(i);
                TxBitsLLR((i-1)*ModType+2) = abs(TxBitsLLR((i-1)*ModType+1)) - 8; % 8 
                TxBitsLLR((i-1)*ModType+3) = abs(TxBitsLLR((i-1)*ModType+2)) - 4;% 4 
                TxBitsLLR((i-1)*ModType+4) = abs(TxBitsLLR((i-1)*ModType+3)) - 2;% 2
    
                TxBitsLLR((i-1)*ModType+5) = symLLR(i+TxAntNum);
                TxBitsLLR((i-1)*ModType+6) = abs(TxBitsLLR((i-1)*ModType+5)) - 8;
                TxBitsLLR((i-1)*ModType+7) = abs(TxBitsLLR((i-1)*ModType+6)) - 4;
                TxBitsLLR((i-1)*ModType+8) = abs(TxBitsLLR((i-1)*ModType+7)) - 2;
            end

    end


end

    


end