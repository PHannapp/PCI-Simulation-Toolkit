import numpy as np
# reference -------------------------------------------------------------------------------
def GHSERZR(T):
    return -7827.595 + 125.64905 * T - 24.1618 * T * np.log(T) - 4.37791E-3 * T**2 + 34971 * T**(-1)
def GHSERTI(T):
    return -8059.921 + 133.615208 * T - 23.9933 * T * np.log(T) - 4.777975E-3 * T**2 + 0.106716E-6 * T**3 + 72636 * T**(-1)
def GHSERH2(T):
    return -9522.97393+78.5273873*T-31.35707*T*np.log(T)+0.0027589925*T**2-0.000000746390667*T**3+56582.3*T**(-1)
def GHSERMN(T):
    return -8115.28 + 130.059*T - 23.4582*T*np.log(T) - 0.00734768*T**2 + 69827*T**(-1)

# Endmembers -------------------------------------------------------------------------------
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
def GTI_MN_H_H_V(T):
   return 1*GHSERTI(T) + 2*GHSERMN(T) + 1.0*GHSERH2(T) + (-62497) + (130.0 * T)
def GTI_MN_H_H_H(T):
   return 1*GHSERTI(T) + 2*GHSERMN(T) + 2.0*GHSERH2(T) + (-75681) + (260.0 * T)
def GTI_MN_H_V_H(T):
   return 1*GHSERTI(T) + 2*GHSERMN(T) + 1.5*GHSERH2(T) + (-77714) + (195.0 * T)
def GTI_MN_H_V_V(T):
   return 1*GHSERTI(T) + 2*GHSERMN(T) + 0.5*GHSERH2(T) + (-67224) + (65.0 * T)
def GZR_MN_V_H_H(T):
   return 1*GHSERZR(T) + 2*GHSERMN(T) + 1.5*GHSERH2(T) + (-103736) + (195.0 * T)
def GZR_MN_V_H_V(T):
   return 1*GHSERZR(T) + 2*GHSERMN(T) + 0.5*GHSERH2(T) + (-66367) + (65.0 * T)
def GZR_MN_V_V_H(T):
   return 1*GHSERZR(T) + 2*GHSERMN(T) + 1.0*GHSERH2(T) + (-76721) + (130.0 * T)
def GZR_MN_V_V_V(T):
   return 1*GHSERZR(T) + 2*GHSERMN(T) + 0*GHSERH2(T) + (-54153) + (0 * T)

# Interaction parameters -------------------------------------------------------------------
def LTIZR_MN_H_H_H_0(T):
   return 41573
def LTIZR_MN_H_H_H_1(T):
   return 43726
def LTIZR_MN_H_H_H_2(T):
   return 28534
def LTIZR_MN_H_H_V_0(T):
   return 26710
def LTIZR_MN_H_H_V_1(T):
   return -21241
def LTIZR_MN_H_H_V_2(T):
   return 47628
def LTIZR_MN_H_V_H_0(T):
   return 31741
def LTIZR_MN_H_V_H_1(T):
   return 3509
def LTIZR_MN_H_V_H_2(T):
   return -80301
def LTIZR_MN_H_V_V_0(T):
   return 16900
def LTIZR_MN_H_V_V_1(T):
   return -10187
def LTIZR_MN_H_V_V_2(T):
   return -4316
def LTIZR_MN_V_H_H_0(T):
   return 30517
def LTIZR_MN_V_H_H_1(T):
   return 18179
def LTIZR_MN_V_H_H_2(T):
   return -47235
def LTIZR_MN_V_H_V_0(T):
   return 10753
def LTIZR_MN_V_H_V_1(T):
   return 5136
def LTIZR_MN_V_H_V_2(T):
   return -13063
def LTIZR_MN_V_V_H_0(T):
   return 26424
def LTIZR_MN_V_V_H_1(T):
   return 18751
def LTIZR_MN_V_V_H_2(T):
   return 13567
def LTIZR_MN_V_V_V_0(T):
   return 12182
def LTIZR_MN_V_V_V_1(T):
   return 4549
def LTIZR_MN_V_V_V_2(T):
   return -2938
