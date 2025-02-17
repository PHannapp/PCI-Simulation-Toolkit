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
def LTI_MN_HV_H_H_0(T):
   return 20799
def LTI_MN_HV_H_H_1(T):
   return -20425
def LTI_MN_HV_H_H_2(T):
   return 61823
def LTI_MN_HV_H_V_0(T):
   return 13494
def LTI_MN_HV_H_V_1(T):
   return 5707
def LTI_MN_HV_H_V_2(T):
   return -21928
def LTI_MN_HV_V_H_0(T):
   return 12692
def LTI_MN_HV_V_H_1(T):
   return -19760
def LTI_MN_HV_V_H_2(T):
   return -47479
def LTI_MN_HV_V_V_0(T):
   return 15628
def LTI_MN_HV_V_V_1(T):
   return -701
def LTI_MN_HV_V_V_2(T):
   return -12008
def LZR_MN_HV_H_H_0(T):
   return -7135
def LZR_MN_HV_H_H_1(T):
   return 31273
def LZR_MN_HV_H_H_2(T):
   return 71151
def LZR_MN_HV_H_V_0(T):
   return -34822
def LZR_MN_HV_H_V_1(T):
   return -24607
def LZR_MN_HV_H_V_2(T):
   return 52217
def LZR_MN_HV_V_H_0(T):
   return -7466
def LZR_MN_HV_V_H_1(T):
   return -1427
def LZR_MN_HV_V_H_2(T):
   return 80325
def LZR_MN_HV_V_V_0(T):
   return 245
def LZR_MN_HV_V_V_1(T):
   return -4083
def LZR_MN_HV_V_V_2(T):
   return 39319
def LTI_MN_H_HV_H_0(T):
   return 16914
def LTI_MN_H_HV_H_1(T):
   return 22931
def LTI_MN_H_HV_H_2(T):
   return 191439
def LTI_MN_H_HV_V_0(T):
   return 38190
def LTI_MN_H_HV_V_1(T):
   return -69799
def LTI_MN_H_HV_V_2(T):
   return -11991
def LTI_MN_V_HV_H_0(T):
   return -8613
def LTI_MN_V_HV_H_1(T):
   return -8933
def LTI_MN_V_HV_H_2(T):
   return 1482
def LTI_MN_V_HV_V_0(T):
   return 17202
def LTI_MN_V_HV_V_1(T):
   return -1102
def LTI_MN_V_HV_V_2(T):
   return 604
def LZR_MN_H_HV_H_0(T):
   return 4832
def LZR_MN_H_HV_H_1(T):
   return 42937
def LZR_MN_H_HV_H_2(T):
   return 3636
def LZR_MN_H_HV_V_0(T):
   return -4975
def LZR_MN_H_HV_V_1(T):
   return -51609
def LZR_MN_H_HV_V_2(T):
   return 4988
def LZR_MN_V_HV_H_0(T):
   return -1478
def LZR_MN_V_HV_H_1(T):
   return -14235
def LZR_MN_V_HV_H_2(T):
   return -122258
def LZR_MN_V_HV_V_0(T):
   return 2841
def LZR_MN_V_HV_V_1(T):
   return -12427
def LZR_MN_V_HV_V_2(T):
   return 9884
def LTI_MN_H_H_HV_0(T):
   return -24391
def LTI_MN_H_H_HV_1(T):
   return -13985
def LTI_MN_H_H_HV_2(T):
   return 71872
def LTI_MN_H_V_HV_0(T):
   return 16964
def LTI_MN_H_V_HV_1(T):
   return -3197
def LTI_MN_H_V_HV_2(T):
   return 46430
def LTI_MN_V_H_HV_0(T):
   return -7473
def LTI_MN_V_H_HV_1(T):
   return -323
def LTI_MN_V_H_HV_2(T):
   return 23892
def LTI_MN_V_V_HV_0(T):
   return 22389
def LTI_MN_V_V_HV_1(T):
   return -19717
def LTI_MN_V_V_HV_2(T):
   return 31760
def LZR_MN_H_H_HV_0(T):
   return -25634
def LZR_MN_H_H_HV_1(T):
   return 39617
def LZR_MN_H_H_HV_2(T):
   return -94264
def LZR_MN_H_V_HV_0(T):
   return -12101
def LZR_MN_H_V_HV_1(T):
   return -14310
def LZR_MN_H_V_HV_2(T):
   return -7503
def LZR_MN_V_H_HV_0(T):
   return -8125
def LZR_MN_V_H_HV_1(T):
   return 4455
def LZR_MN_V_H_HV_2(T):
   return -20746
def LZR_MN_V_V_HV_0(T):
   return 1399
def LZR_MN_V_V_HV_1(T):
   return -6853
def LZR_MN_V_V_HV_2(T):
   return 16911
