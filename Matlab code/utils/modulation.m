function [TxSymbol] = modulation(TxBits,ModType,s)
%   Detailed explanation goes here

half_sym=ModType/2;

TxBits1=reshape(TxBits,ModType,[]);
TxBit_real=TxBits1(1:ModType/2,:);  % real part modulation
TxBit_imag=TxBits1(ModType/2+1:end,:); % imag part modulation
TxBit_real_int=bit2int(TxBit_real,half_sym)+1;
TxBit_imag_int=bit2int(TxBit_imag,half_sym)+1;

TxSymbol=s(TxBit_real_int)+1i*s(TxBit_imag_int);

end