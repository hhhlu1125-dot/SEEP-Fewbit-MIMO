
%% Set MIMO Case  32*8 16QAM 2-bit;
clear;
clc;

N_t = 32;
N_r = 128;
A = 8;
I_max = 5;
L_net = 10;
L_QVB = 10;

k = 14;
m = 8;

EP_NRVs = N_r*N_t^2 + 2/3*N_t^3 + I_max*(2/3*N_t^3+4*N_t*A);
% EP_Add = N_r*N_t^2 + 2/3*N_t^3 + I_max*(2/3*N_t^3+5*N_t*A);
% EP_Exp = I_max*N_t*A;

BEP_NRVs = N_r^2*N_t + N_r*N_t^2 + 2/3*N_t^3 + I_max*(2/3*N_t^3+4*N_t*A);
% BEP_Add = N_r^2*N_t + N_r*N_t^2 + 2/3*N_t^3 + I_max*(2/3*N_t^3+5*N_t*A);
% BEP_Exp = I_max*N_t*A;

SEEP_NRVs = N_r^2*N_t + 2*N_r*N_t^2 + 2/3*N_t^3 + I_max*(N_r^2*N_t+2*N_r*N_t^2+4/3*N_t^3+4*N_t*A);
% SEEP_Add = N_r^2*N_t + 2*N_r*N_t^2 + 2/3*N_t^3 + I_max*(N_r^2*N_t+2*N_r*N_t^2+4/3*N_t^3+5*N_t*A);
% SEEP_Exp = I_max*(N_t*A+4*N_r);

aSEEP_NRVs = 2*N_r*N_t^2 + 2/3*N_t^3 + I_max*(2*N_r*N_t^2+2/3*N_t^3+4*N_t*A);
% aSEEP_Add = 2*N_r*N_t^2 + 2/3*N_t^3 + I_max*(2*N_r*N_t^2+2/3*N_t^3+5*N_t*A);
% aSEEP_Exp = I_max*(N_t*A+4*N_r);

% CFxPOs_EP = k^2*EP_Mul+k*EP_Add+m*k^2*EP_Exp;
% CFxPOs_BEP = k^2*BEP_Mul+k*BEP_Add+m*k^2*BEP_Exp;
% CFxPOs_SEEP = k^2*SEEP_Mul+k*SEEP_Add+m*k^2*SEEP_Exp;
% CFxPOs_aSEEP = k^2*aSEEP_Mul+k*aSEEP_Add+m*k^2*aSEEP_Exp;
% CReduction_EP = abs((CFxPOs_EP-CFxPOs_aSEEP))/CFxPOs_EP;
% CReduction_BEP = abs((CFxPOs_BEP-CFxPOs_aSEEP))/CFxPOs_BEP;
% CReduction_SEEP = abs((CFxPOs_SEEP-CFxPOs_aSEEP))/CFxPOs_SEEP;

a = SEEP_NRVs / EP_NRVs;
b = SEEP_NRVs / BEP_NRVs;
c = (SEEP_NRVs-aSEEP_NRVs)/ SEEP_NRVs;

FBM_DetNet_NRVs = 2*L_net*N_r*N_t;
LMMSE_QVB_NRVs = L_QVB*(2/3*N_r^3+2*N_r^2*N_t);
