import numpy as np
# reference -------------------------------------------------------------------------------
def GHSERAL(T):
    return -7976.15 + 137.093038*T - 24.3671976*T*np.log(T) - 1.884662E-3*T**2 - 0.877664E-6*T**3 + 74092*T**(-1)
def GHSERCE(T):
    return -7160.519+84.23022*T-22.3664*T*np.log(T)-0.0067103*T**2-0.000000320773*T**3-18117*T**(-1)
def GHSERLA(T):
    return -7968.403 + 120.284604*T - 26.34*T*np.log(T) - 0.001295165*T**2
def GHSERNI(T):
    return -5179.159+117.854*T-22.096*T*np.log(T)-0.0048407*T**2
def GHSERH2(T):
    return -9522.97393+78.5273873*T-31.35707*T*np.log(T)+0.0027589925*T**2-0.000000746390667*T**3+56582.3*T**(-1)

# Endmembers -------------------------------------------------------------------------------
def GCE_AL_H_H(T):
   return 1*GHSERCE(T) + 5*GHSERAL(T) + 3.5*GHSERH2(T) + (310249) + (427.0 * T)
def GCE_AL_H_V(T):
   return 1*GHSERCE(T) + 5*GHSERAL(T) + 0.5*GHSERH2(T) + (-62751) + (61.0 * T)
def GCE_AL_V_H(T):
   return 1*GHSERCE(T) + 5*GHSERAL(T) + 3.0*GHSERH2(T) + (1099926) + (366.0 * T)
def GCE_AL_V_V(T):
   return 1*GHSERCE(T) + 5*GHSERAL(T) + 0*GHSERH2(T) + (113908) + (0 * T)
def GCE_NI_H_H(T):
   return 1*GHSERCE(T) + 5*GHSERNI(T) + 3.5*GHSERH2(T) + (-446010) + (427.0 * T)
def GCE_NI_H_V(T):
   return 1*GHSERCE(T) + 5*GHSERNI(T) + 0.5*GHSERH2(T) + (-83391) + (61.0 * T)
def GCE_NI_V_H(T):
   return 1*GHSERCE(T) + 5*GHSERNI(T) + 3.0*GHSERH2(T) + (310727) + (366.0 * T)
def GCE_NI_V_V(T):
   return 1*GHSERCE(T) + 5*GHSERNI(T) + 0*GHSERH2(T) + (-212469) + (0 * T)
def GLA_AL_H_H(T):
   return 1*GHSERLA(T) + 5*GHSERAL(T) + 3.5*GHSERH2(T) + (338322) + (427.0 * T)
def GLA_AL_H_V(T):
   return 1*GHSERLA(T) + 5*GHSERAL(T) + 0.5*GHSERH2(T) + (-57115) + (61.0 * T)
def GLA_AL_V_H(T):
   return 1*GHSERLA(T) + 5*GHSERAL(T) + 3.0*GHSERH2(T) + (1089440) + (366.0 * T)
def GLA_AL_V_V(T):
   return 1*GHSERLA(T) + 5*GHSERAL(T) + 0*GHSERH2(T) + (187950) + (0 * T)
def GLA_NI_H_H(T):
   return 1*GHSERLA(T) + 5*GHSERNI(T) + 3.5*GHSERH2(T) + (-453534) + (427.0 * T)
def GLA_NI_H_V(T):
   return 1*GHSERLA(T) + 5*GHSERNI(T) + 0.5*GHSERH2(T) + (-116041) + (61.0 * T)
def GLA_NI_V_H(T):
   return 1*GHSERLA(T) + 5*GHSERNI(T) + 3.0*GHSERH2(T) + (292668) + (366.0 * T)
def GLA_NI_V_V(T):
   return 1*GHSERLA(T) + 5*GHSERNI(T) + 0*GHSERH2(T) + (-194872) + (0 * T)

# Interaction parameters -------------------------------------------------------------------
def LCELA_AL_H_H_0(T):
   return -3257
def LCELA_AL_H_H_1(T):
   return 427
def LCELA_AL_H_H_2(T):
   return 346
def LCELA_AL_H_V_0(T):
   return 235382
def LCELA_AL_H_V_1(T):
   return 8060
def LCELA_AL_H_V_2(T):
   return 556142
def LCELA_AL_V_H_0(T):
   return -2640743
def LCELA_AL_V_H_1(T):
   return -100996
def LCELA_AL_V_H_2(T):
   return -6068666
def LCELA_AL_V_V_0(T):
   return 18
def LCELA_AL_V_V_1(T):
   return 186
def LCELA_AL_V_V_2(T):
   return -485
def LCELA_NI_H_H_0(T):
   return -2731
def LCELA_NI_H_H_1(T):
   return -109
def LCELA_NI_H_H_2(T):
   return -407
def LCELA_NI_H_V_0(T):
   return -3649
def LCELA_NI_H_V_1(T):
   return 1210
def LCELA_NI_H_V_2(T):
   return -7222
def LCELA_NI_V_H_0(T):
   return -1996396
def LCELA_NI_V_H_1(T):
   return -17094
def LCELA_NI_V_H_2(T):
   return -4521711
def LCELA_NI_V_V_0(T):
   return -1758
def LCELA_NI_V_V_1(T):
   return 36
def LCELA_NI_V_V_2(T):
   return -1122
def LCE_ALNI_H_H_0(T):
   return -783815
def LCE_ALNI_H_H_1(T):
   return 86400
def LCE_ALNI_H_H_2(T):
   return 101856
def LCE_ALNI_H_V_0(T):
   return -686278
def LCE_ALNI_H_V_1(T):
   return 51013
def LCE_ALNI_H_V_2(T):
   return 235874
def LCE_ALNI_V_H_0(T):
   return -1409783
def LCE_ALNI_V_H_1(T):
   return 181223
def LCE_ALNI_V_H_2(T):
   return -3707248
def LCE_ALNI_V_V_0(T):
   return -879060
def LCE_ALNI_V_V_1(T):
   return -105421
def LCE_ALNI_V_V_2(T):
   return 92771
def LLA_ALNI_H_H_0(T):
   return -839331
def LLA_ALNI_H_H_1(T):
   return 50272
def LLA_ALNI_H_H_2(T):
   return 94326
def LLA_ALNI_H_V_0(T):
   return -717837
def LLA_ALNI_H_V_1(T):
   return 95333
def LLA_ALNI_H_V_2(T):
   return 307221
def LLA_ALNI_V_H_0(T):
   return -1131782
def LLA_ALNI_V_H_1(T):
   return -319226
def LLA_ALNI_V_H_2(T):
   return -3165085
def LLA_ALNI_V_V_0(T):
   return -868128
def LLA_ALNI_V_V_1(T):
   return -102306
def LLA_ALNI_V_V_2(T):
   return 82009
def LCE_AL_HV_H_0(T):
   return -1873856
def LCE_AL_HV_H_1(T):
   return 3317954
def LCE_AL_HV_H_2(T):
   return -4562136
def LCE_AL_HV_V_0(T):
   return -503983
def LCE_AL_HV_V_1(T):
   return 1327801
def LCE_AL_HV_V_2(T):
   return -1140533
def LCE_NI_HV_H_0(T):
   return -1453443
def LCE_NI_HV_H_1(T):
   return 2533922
def LCE_NI_HV_H_2(T):
   return -3187986
def LCE_NI_HV_V_0(T):
   return 413716
def LCE_NI_HV_V_1(T):
   return -745933
def LCE_NI_HV_V_2(T):
   return 981383
def LLA_AL_HV_H_0(T):
   return -1828733
def LLA_AL_HV_H_1(T):
   return 3256002
def LLA_AL_HV_H_2(T):
   return -4034031
def LLA_AL_HV_V_0(T):
   return -623867
def LLA_AL_HV_V_1(T):
   return 1543551
def LLA_AL_HV_V_2(T):
   return -1608962
def LLA_NI_HV_H_0(T):
   return -1429032
def LLA_NI_HV_H_1(T):
   return 2515557
def LLA_NI_HV_H_2(T):
   return -3247603
def LLA_NI_HV_V_0(T):
   return 297481
def LLA_NI_HV_V_1(T):
   return -551999
def LLA_NI_HV_V_2(T):
   return 715982
def LCE_AL_H_HV_0(T):
   return 448590
def LCE_AL_H_HV_1(T):
   return -598690
def LCE_AL_H_HV_2(T):
   return 327393
def LCE_AL_V_HV_0(T):
   return -2404063
def LCE_AL_V_HV_1(T):
   return -1453716
def LCE_AL_V_HV_2(T):
   return -6658202
def LCE_NI_H_HV_0(T):
   return -225162
def LCE_NI_H_HV_1(T):
   return 66439
def LCE_NI_H_HV_2(T):
   return -61566
def LCE_NI_V_HV_0(T):
   return -895872
def LCE_NI_V_HV_1(T):
   return -3490651
def LCE_NI_V_HV_2(T):
   return -1401334
def LLA_AL_H_HV_0(T):
   return 395698
def LLA_AL_H_HV_1(T):
   return -587029
def LLA_AL_H_HV_2(T):
   return 860080
def LLA_AL_V_HV_0(T):
   return -2497648
def LLA_AL_V_HV_1(T):
   return -660720
def LLA_AL_V_HV_2(T):
   return -7198891
def LLA_NI_H_HV_0(T):
   return -244432
def LLA_NI_H_HV_1(T):
   return 83670
def LLA_NI_H_HV_2(T):
   return -28642
def LLA_NI_V_HV_0(T):
   return -1096629
def LLA_NI_V_HV_1(T):
   return -3124236
def LLA_NI_V_HV_2(T):
   return -1711360
