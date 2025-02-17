import numpy as np
# reference -------------------------------------------------------------------------------
def GHSERTI(T):
    return -8059.921 + 133.615208 * T - 23.9933 * T * np.log(T) - 4.777975E-3 * T**2 + 0.106716E-6 * T**3 + 72636 * T**(-1)
def GHSERH2(T):
    return -9522.97393+78.5273873*T-31.35707*T*np.log(T)+0.0027589925*T**2-0.000000746390667*T**3+56582.3*T**(-1)
def GHSERZR(T):
    return -7827.595 + 125.64905 * T - 24.1618 * T * np.log(T) - 4.37791E-3 * T**2 + 34971 * T**(-1)
def GHSERMN(T):
    return -8115.28 + 130.059*T - 23.4582*T*np.log(T) - 0.00734768*T**2 + 69827*T**(-1)

# Endmembers -------------------------------------------------------------------------------
def GZR_MN_V_V_H(T):
   return 1*GHSERZR(T) + 2*GHSERMN(T) + 1.0*GHSERH2(T) + (-76721) + (130.0 * T)
def GZR_MN_V_V_V(T):
   return 1*GHSERZR(T) + 2*GHSERMN(T) + 0*GHSERH2(T) + (-54153) + (0 * T)
def GTI_MN_H_V_H(T):
   return 1*GHSERTI(T) + 2*GHSERMN(T) + 1.5*GHSERH2(T) + (-77714) + (195.0 * T)
def GTI_MN_H_V_V(T):
   return 1*GHSERTI(T) + 2*GHSERMN(T) + 0.5*GHSERH2(T) + (-67224) + (65.0 * T)
def GZR_MN_V_H_H(T):
   return 1*GHSERZR(T) + 2*GHSERMN(T) + 1.5*GHSERH2(T) + (-103736) + (195.0 * T)
def GZR_MN_V_H_V(T):
   return 1*GHSERZR(T) + 2*GHSERMN(T) + 0.5*GHSERH2(T) + (-66367) + (65.0 * T)
def GTI_MN_H_H_V(T):
   return 1*GHSERTI(T) + 2*GHSERMN(T) + 1.0*GHSERH2(T) + (-62497) + (130.0 * T)
def GTI_MN_H_H_H(T):
   return 1*GHSERTI(T) + 2*GHSERMN(T) + 2.0*GHSERH2(T) + (-75681) + (260.0 * T)
def GTI_MN_V_H_H(T):
   return 1*GHSERTI(T) + 2*GHSERMN(T) + 1.5*GHSERH2(T) + (-66043) + (195.0 * T)
def GTI_MN_V_H_V(T):
   return 1*GHSERTI(T) + 2*GHSERMN(T) + 0.5*GHSERH2(T) + (-62165) + (65.0 * T)
def GTI_MN_V_V_H(T):
   return 1*GHSERTI(T) + 2*GHSERMN(T) + 1.0*GHSERH2(T) + (-62343) + (130.0 * T)
def GTI_MN_V_V_V(T):
   return 1*GHSERTI(T) + 2*GHSERMN(T) + 0*GHSERH2(T) + (-76200) + (0 * T)
def GZR_MN_H_H_H(T):
   return 1*GHSERZR(T) + 2*GHSERMN(T) + 2.0*GHSERH2(T) + (-123674) + (260.0 * T)
def GZR_MN_H_H_V(T):
   return 1*GHSERZR(T) + 2*GHSERMN(T) + 1.0*GHSERH2(T) + (-72711) + (130.0 * T)
def GZR_MN_H_V_H(T):
   return 1*GHSERZR(T) + 2*GHSERMN(T) + 1.5*GHSERH2(T) + (-97596) + (195.0 * T)
def GZR_MN_H_V_V(T):
   return 1*GHSERZR(T) + 2*GHSERMN(T) + 0.5*GHSERH2(T) + (-63872) + (65.0 * T)

# Interaction parameters -------------------------------------------------------------------
def LTIZR_MN_H_H_H_0(T):
   return 45157
def LTIZR_MN_H_H_H_1(T):
   return 0
def LTIZR_MN_H_H_H_2(T):
   return 0
def LTIZR_MN_H_H_V_0(T):
   return 33014
def LTIZR_MN_H_H_V_1(T):
   return 0
def LTIZR_MN_H_H_V_2(T):
   return 0
def LTIZR_MN_H_V_H_0(T):
   return 21028
def LTIZR_MN_H_V_H_1(T):
   return 0
def LTIZR_MN_H_V_H_2(T):
   return 0
def LTIZR_MN_H_V_V_0(T):
   return 16336
def LTIZR_MN_H_V_V_1(T):
   return 0
def LTIZR_MN_H_V_V_2(T):
   return 0
def LTIZR_MN_V_H_H_0(T):
   return 24269
def LTIZR_MN_V_H_H_1(T):
   return 0
def LTIZR_MN_V_H_H_2(T):
   return 0
def LTIZR_MN_V_H_V_0(T):
   return 9024
def LTIZR_MN_V_H_V_1(T):
   return 0
def LTIZR_MN_V_H_V_2(T):
   return 0
def LTIZR_MN_V_V_H_0(T):
   return 28220
def LTIZR_MN_V_V_H_1(T):
   return 0
def LTIZR_MN_V_V_H_2(T):
   return 0
def LTIZR_MN_V_V_V_0(T):
   return 11793
def LTIZR_MN_V_V_V_1(T):
   return 0
def LTIZR_MN_V_V_V_2(T):
   return 0
def LTI_MN_HV_H_H_0(T):
   return 28872
def LTI_MN_HV_H_H_1(T):
   return 0
def LTI_MN_HV_H_H_2(T):
   return 0
def LTI_MN_HV_H_V_0(T):
   return 10592
def LTI_MN_HV_H_V_1(T):
   return 0
def LTI_MN_HV_H_V_2(T):
   return 0
def LTI_MN_HV_V_H_0(T):
   return 6365
def LTI_MN_HV_V_H_1(T):
   return 0
def LTI_MN_HV_V_H_2(T):
   return 0
def LTI_MN_HV_V_V_0(T):
   return 14039
def LTI_MN_HV_V_V_1(T):
   return 0
def LTI_MN_HV_V_V_2(T):
   return 0
def LZR_MN_HV_H_H_0(T):
   return 2505
def LZR_MN_HV_H_H_1(T):
   return 0
def LZR_MN_HV_H_H_2(T):
   return 0
def LZR_MN_HV_H_V_0(T):
   return -27934
def LZR_MN_HV_H_V_1(T):
   return 0
def LZR_MN_HV_H_V_2(T):
   return 0
def LZR_MN_HV_V_H_0(T):
   return 3163
def LZR_MN_HV_V_H_1(T):
   return 0
def LZR_MN_HV_V_H_2(T):
   return 0
def LZR_MN_HV_V_V_0(T):
   return 5449
def LZR_MN_HV_V_V_1(T):
   return 0
def LZR_MN_HV_V_V_2(T):
   return 0
def LTI_MN_H_HV_H_0(T):
   return 41931
def LTI_MN_H_HV_H_1(T):
   return 0
def LTI_MN_H_HV_H_2(T):
   return 0
def LTI_MN_H_HV_V_0(T):
   return 36963
def LTI_MN_H_HV_V_1(T):
   return 0
def LTI_MN_H_HV_V_2(T):
   return 0
def LTI_MN_V_HV_H_0(T):
   return -8412
def LTI_MN_V_HV_H_1(T):
   return 0
def LTI_MN_V_HV_H_2(T):
   return 0
def LTI_MN_V_HV_V_0(T):
   return 17282
def LTI_MN_V_HV_V_1(T):
   return 0
def LTI_MN_V_HV_V_2(T):
   return 0
def LZR_MN_H_HV_H_0(T):
   return 5391
def LZR_MN_H_HV_H_1(T):
   return 0
def LZR_MN_H_HV_H_2(T):
   return 0
def LZR_MN_H_HV_V_0(T):
   return -4467
def LZR_MN_H_HV_V_1(T):
   return 0
def LZR_MN_H_HV_V_2(T):
   return 0
def LZR_MN_V_HV_H_0(T):
   return -17747
def LZR_MN_V_HV_H_1(T):
   return 0
def LZR_MN_V_HV_H_2(T):
   return 0
def LZR_MN_V_HV_V_0(T):
   return 4150
def LZR_MN_V_HV_V_1(T):
   return 0
def LZR_MN_V_HV_V_2(T):
   return 0
def LTI_MN_H_H_HV_0(T):
   return -15147
def LTI_MN_H_H_HV_1(T):
   return 0
def LTI_MN_H_H_HV_2(T):
   return 0
def LTI_MN_H_V_HV_0(T):
   return 23707
def LTI_MN_H_V_HV_1(T):
   return 0
def LTI_MN_H_V_HV_2(T):
   return 0
def LTI_MN_V_H_HV_0(T):
   return -4076
def LTI_MN_V_H_HV_1(T):
   return 0
def LTI_MN_V_H_HV_2(T):
   return 0
def LTI_MN_V_V_HV_0(T):
   return 27008
def LTI_MN_V_V_HV_1(T):
   return 0
def LTI_MN_V_V_HV_2(T):
   return 0
def LZR_MN_H_H_HV_0(T):
   return -39001
def LZR_MN_H_H_HV_1(T):
   return 0
def LZR_MN_H_H_HV_2(T):
   return 0
def LZR_MN_H_V_HV_0(T):
   return -13090
def LZR_MN_H_V_HV_1(T):
   return 0
def LZR_MN_H_V_HV_2(T):
   return 0
def LZR_MN_V_H_HV_0(T):
   return -11076
def LZR_MN_V_H_HV_1(T):
   return 0
def LZR_MN_V_H_HV_2(T):
   return 0
def LZR_MN_V_V_HV_0(T):
   return 3811
def LZR_MN_V_V_HV_1(T):
   return 0
def LZR_MN_V_V_HV_2(T):
   return 0
