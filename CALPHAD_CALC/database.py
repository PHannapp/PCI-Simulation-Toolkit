import numpy as np
# reference -------------------------------------------------------------------------------
def GHSERLA(T):
    return -7968.403 + 120.284604*T - 26.34*T*np.log(T) - 0.001295165*T**2
def GHSERNI(T):
    return -5179.159+117.854*T-22.096*T*np.log(T)-0.0048407*T**2
def GHSERCE(T):
    return -7160.519+84.23022*T-22.3664*T*np.log(T)-0.0067103*T**2-0.000000320773*T**3-18117*T**(-1)
def GHSERH2(T):
    return -9522.97393+78.5273873*T-31.35707*T*np.log(T)+0.0027589925*T**2-0.000000746390667*T**3+56582.3*T**(-1)
def GHSERAL(T):
    return -7976.15 + 137.093038*T - 24.3671976*T*np.log(T) - 1.884662E-3*T**2 - 0.877664E-6*T**3 + 74092*T**(-1)

# Endmembers -------------------------------------------------------------------------------
def GLA_NI_NI_H_H(T):
   return 1*GHSERLA(T) + 2*GHSERNI(T) + 3*GHSERNI(T) + 3.5*GHSERH2(T) + (-320282) + (455.0 * T)
def GLA_NI_NI_H_V(T):
   return 1*GHSERLA(T) + 2*GHSERNI(T) + 3*GHSERNI(T) + 0.5*GHSERH2(T) + (-169929) + (65.0 * T)
def GCE_NI_NI_H_H(T):
   return 1*GHSERCE(T) + 2*GHSERNI(T) + 3*GHSERNI(T) + 3.5*GHSERH2(T) + (-295649) + (455.0 * T)
def GCE_NI_NI_H_V(T):
   return 1*GHSERCE(T) + 2*GHSERNI(T) + 3*GHSERNI(T) + 0.5*GHSERH2(T) + (-185493) + (65.0 * T)
def GLA_NI_NI_V_H(T):
   return 1*GHSERLA(T) + 2*GHSERNI(T) + 3*GHSERNI(T) + 3.0*GHSERH2(T) + (-270252) + (390.0 * T)
def GLA_NI_NI_V_V(T):
   return 1*GHSERLA(T) + 2*GHSERNI(T) + 3*GHSERNI(T) + 0*GHSERH2(T) + (-177701) + (0 * T)
def GCE_AL_AL_V_H(T):
   return 1*GHSERCE(T) + 2*GHSERAL(T) + 3*GHSERAL(T) + 3.0*GHSERH2(T) + (341049) + (390.0 * T)
def GCE_AL_AL_V_V(T):
   return 1*GHSERCE(T) + 2*GHSERAL(T) + 3*GHSERAL(T) + 0*GHSERH2(T) + (-108091) + (0 * T)
def GLA_AL_NI_H_H(T):
   return 1*GHSERLA(T) + 2*GHSERAL(T) + 3*GHSERNI(T) + 3.5*GHSERH2(T) + (57619) + (455.0 * T)
def GLA_AL_NI_H_V(T):
   return 1*GHSERLA(T) + 2*GHSERAL(T) + 3*GHSERNI(T) + 0.5*GHSERH2(T) + (-212054) + (65.0 * T)
def GLA_AL_NI_V_H(T):
   return 1*GHSERLA(T) + 2*GHSERAL(T) + 3*GHSERNI(T) + 3.0*GHSERH2(T) + (-56922) + (390.0 * T)
def GLA_AL_NI_V_V(T):
   return 1*GHSERLA(T) + 2*GHSERAL(T) + 3*GHSERNI(T) + 0*GHSERH2(T) + (-319839) + (0 * T)
def GLA_NI_AL_V_H(T):
   return 1*GHSERLA(T) + 2*GHSERNI(T) + 3*GHSERAL(T) + 3.0*GHSERH2(T) + (-133337) + (390.0 * T)
def GLA_NI_AL_V_V(T):
   return 1*GHSERLA(T) + 2*GHSERNI(T) + 3*GHSERAL(T) + 0*GHSERH2(T) + (-282404) + (0 * T)
def GCE_NI_AL_H_V(T):
   return 1*GHSERCE(T) + 2*GHSERNI(T) + 3*GHSERAL(T) + 0.5*GHSERH2(T) + (-195327) + (65.0 * T)
def GLA_NI_AL_H_H(T):
   return 1*GHSERLA(T) + 2*GHSERNI(T) + 3*GHSERAL(T) + 3.5*GHSERH2(T) + (-66799) + (455.0 * T)
def GLA_NI_AL_H_V(T):
   return 1*GHSERLA(T) + 2*GHSERNI(T) + 3*GHSERAL(T) + 0.5*GHSERH2(T) + (-217033) + (65.0 * T)
def GCE_NI_AL_H_H(T):
   return 1*GHSERCE(T) + 2*GHSERNI(T) + 3*GHSERAL(T) + 3.5*GHSERH2(T) + (-25512) + (455.0 * T)
def GLA_AL_AL_V_H(T):
   return 1*GHSERLA(T) + 2*GHSERAL(T) + 3*GHSERAL(T) + 3.0*GHSERH2(T) + (282004) + (390.0 * T)
def GLA_AL_AL_V_V(T):
   return 1*GHSERLA(T) + 2*GHSERAL(T) + 3*GHSERAL(T) + 0*GHSERH2(T) + (-126419) + (0 * T)
def GCE_NI_NI_V_H(T):
   return 1*GHSERCE(T) + 2*GHSERNI(T) + 3*GHSERNI(T) + 3.0*GHSERH2(T) + (-259566) + (390.0 * T)
def GCE_NI_NI_V_V(T):
   return 1*GHSERCE(T) + 2*GHSERNI(T) + 3*GHSERNI(T) + 0*GHSERH2(T) + (-216483) + (0 * T)
def GCE_AL_NI_V_H(T):
   return 1*GHSERCE(T) + 2*GHSERAL(T) + 3*GHSERNI(T) + 3.0*GHSERH2(T) + (-27057) + (390.0 * T)
def GCE_AL_NI_V_V(T):
   return 1*GHSERCE(T) + 2*GHSERAL(T) + 3*GHSERNI(T) + 0*GHSERH2(T) + (-319839) + (0 * T)
def GCE_NI_AL_V_H(T):
   return 1*GHSERCE(T) + 2*GHSERNI(T) + 3*GHSERAL(T) + 3.0*GHSERH2(T) + (-88071) + (390.0 * T)
def GCE_NI_AL_V_V(T):
   return 1*GHSERCE(T) + 2*GHSERNI(T) + 3*GHSERAL(T) + 0*GHSERH2(T) + (-276139) + (0 * T)
def GCE_AL_AL_H_V(T):
   return 1*GHSERCE(T) + 2*GHSERAL(T) + 3*GHSERAL(T) + 0.5*GHSERH2(T) + (21183) + (65.0 * T)
def GCE_AL_AL_H_H(T):
   return 1*GHSERCE(T) + 2*GHSERAL(T) + 3*GHSERAL(T) + 3.5*GHSERH2(T) + (432302) + (455.0 * T)
def GLA_AL_AL_H_H(T):
   return 1*GHSERLA(T) + 2*GHSERAL(T) + 3*GHSERAL(T) + 3.5*GHSERH2(T) + (142340) + (455.0 * T)
def GLA_AL_AL_H_V(T):
   return 1*GHSERLA(T) + 2*GHSERAL(T) + 3*GHSERAL(T) + 0.5*GHSERH2(T) + (-10527) + (65.0 * T)
def GCE_AL_NI_H_V(T):
   return 1*GHSERCE(T) + 2*GHSERAL(T) + 3*GHSERNI(T) + 0.5*GHSERH2(T) + (-196670) + (65.0 * T)
def GCE_AL_NI_H_H(T):
   return 1*GHSERCE(T) + 2*GHSERAL(T) + 3*GHSERNI(T) + 3.5*GHSERH2(T) + (94350) + (455.0 * T)

# Interaction parameters -------------------------------------------------------------------
def LCELA_AL_AL_H_H_0(T):
   return -762468
def LCELA_AL_AL_H_H_1(T):
   return -1744568
def LCELA_AL_AL_H_H_2(T):
   return 0
def LCELA_AL_AL_H_V_0(T):
   return -74797
def LCELA_AL_AL_H_V_1(T):
   return 56242
def LCELA_AL_AL_H_V_2(T):
   return 0
def LCELA_AL_AL_V_H_0(T):
   return -1434059
def LCELA_AL_AL_V_H_1(T):
   return -264395
def LCELA_AL_AL_V_H_2(T):
   return 0
def LCELA_AL_AL_V_V_0(T):
   return 72929
def LCELA_AL_AL_V_V_1(T):
   return 377275
def LCELA_AL_AL_V_V_2(T):
   return 0
def LCELA_AL_NI_H_H_0(T):
   return -1480206
def LCELA_AL_NI_H_H_1(T):
   return -370642
def LCELA_AL_NI_H_H_2(T):
   return 0
def LCELA_AL_NI_H_V_0(T):
   return -285454
def LCELA_AL_NI_H_V_1(T):
   return 1050530
def LCELA_AL_NI_H_V_2(T):
   return 0
def LCELA_AL_NI_V_H_0(T):
   return -967350
def LCELA_AL_NI_V_H_1(T):
   return -255573
def LCELA_AL_NI_V_H_2(T):
   return 0
def LCELA_AL_NI_V_V_0(T):
   return -60151
def LCELA_AL_NI_V_V_1(T):
   return 386806
def LCELA_AL_NI_V_V_2(T):
   return 0
def LCELA_NI_AL_H_H_0(T):
   return -504119
def LCELA_NI_AL_H_H_1(T):
   return -1611524
def LCELA_NI_AL_H_H_2(T):
   return 0
def LCELA_NI_AL_H_V_0(T):
   return -449732
def LCELA_NI_AL_H_V_1(T):
   return -22027
def LCELA_NI_AL_H_V_2(T):
   return 0
def LCELA_NI_AL_V_H_0(T):
   return -895374
def LCELA_NI_AL_V_H_1(T):
   return -326932
def LCELA_NI_AL_V_H_2(T):
   return 0
def LCELA_NI_AL_V_V_0(T):
   return 1699651
def LCELA_NI_AL_V_V_1(T):
   return -3006584
def LCELA_NI_AL_V_V_2(T):
   return 0
def LCELA_NI_NI_H_H_0(T):
   return -21835
def LCELA_NI_NI_H_H_1(T):
   return -16660
def LCELA_NI_NI_H_H_2(T):
   return 0
def LCELA_NI_NI_H_V_0(T):
   return -150269
def LCELA_NI_NI_H_V_1(T):
   return 75859
def LCELA_NI_NI_H_V_2(T):
   return 0
def LCELA_NI_NI_V_H_0(T):
   return -95953
def LCELA_NI_NI_V_H_1(T):
   return 15739
def LCELA_NI_NI_V_H_2(T):
   return 0
def LCE_ALNI_AL_H_H_0(T):
   return -1658452
def LCE_ALNI_AL_H_H_1(T):
   return -274829
def LCE_ALNI_AL_H_H_2(T):
   return 0
def LCE_ALNI_AL_H_V_0(T):
   return -743087
def LCE_ALNI_AL_H_V_1(T):
   return -288368
def LCE_ALNI_AL_H_V_2(T):
   return 0
def LCE_ALNI_AL_V_H_0(T):
   return -1449895
def LCE_ALNI_AL_V_H_1(T):
   return -392430
def LCE_ALNI_AL_V_H_2(T):
   return 0
def LCE_ALNI_AL_V_V_0(T):
   return -364494
def LCE_ALNI_AL_V_V_1(T):
   return 176537
def LCE_ALNI_AL_V_V_2(T):
   return 0
def LCE_ALNI_NI_H_H_0(T):
   return -1241089
def LCE_ALNI_NI_H_H_1(T):
   return -1079348
def LCE_ALNI_NI_H_H_2(T):
   return 0
def LCE_ALNI_NI_H_V_0(T):
   return -592779
def LCE_ALNI_NI_H_V_1(T):
   return -396580
def LCE_ALNI_NI_H_V_2(T):
   return 0
def LCE_ALNI_NI_V_H_0(T):
   return -1063137
def LCE_ALNI_NI_V_H_1(T):
   return -640762
def LCE_ALNI_NI_V_H_2(T):
   return 0
def LCE_ALNI_NI_V_V_0(T):
   return 4000000
def LCE_ALNI_NI_V_V_1(T):
   return 4000000
def LCE_ALNI_NI_V_V_2(T):
   return 0
def LLA_ALNI_AL_H_H_0(T):
   return -976081
def LLA_ALNI_AL_H_H_1(T):
   return 564351
def LLA_ALNI_AL_H_H_2(T):
   return 0
def LLA_ALNI_AL_H_V_0(T):
   return -759998
def LLA_ALNI_AL_H_V_1(T):
   return 44167
def LLA_ALNI_AL_H_V_2(T):
   return 0
def LLA_ALNI_AL_V_H_0(T):
   return -208700
def LLA_ALNI_AL_V_H_1(T):
   return 2056655
def LLA_ALNI_AL_V_H_2(T):
   return 0
def LLA_ALNI_NI_V_H_0(T):
   return -933263
def LLA_ALNI_NI_V_H_1(T):
   return -3546517
def LLA_ALNI_NI_V_H_2(T):
   return 0
def LCE_AL_ALNI_H_H_0(T):
   return -1880772
def LCE_AL_ALNI_H_H_1(T):
   return 459892
def LCE_AL_ALNI_H_H_2(T):
   return 0
def LCE_AL_ALNI_H_V_0(T):
   return -626409
def LCE_AL_ALNI_H_V_1(T):
   return 42692
def LCE_AL_ALNI_H_V_2(T):
   return 0
def LCE_AL_ALNI_V_H_0(T):
   return -1339748
def LCE_AL_ALNI_V_H_1(T):
   return 18576
def LCE_AL_ALNI_V_H_2(T):
   return 0
def LCE_NI_ALNI_H_H_0(T):
   return -636909
def LCE_NI_ALNI_H_H_1(T):
   return -547661
def LCE_NI_ALNI_H_H_2(T):
   return 0
def LCE_NI_ALNI_H_V_0(T):
   return -2949768
def LCE_NI_ALNI_H_V_1(T):
   return 2926722
def LCE_NI_ALNI_H_V_2(T):
   return 0
def LCE_NI_ALNI_V_H_0(T):
   return -1280778
def LCE_NI_ALNI_V_H_1(T):
   return 2574391
def LCE_NI_ALNI_V_H_2(T):
   return 0
def LCE_NI_ALNI_V_V_0(T):
   return -1814352
def LCE_NI_ALNI_V_V_1(T):
   return -2190661
def LCE_NI_ALNI_V_V_2(T):
   return 0
def LLA_AL_ALNI_H_H_0(T):
   return -1001296
def LLA_AL_ALNI_H_H_1(T):
   return 481074
def LLA_AL_ALNI_H_H_2(T):
   return 0
def LLA_AL_ALNI_H_V_0(T):
   return -693307
def LLA_AL_ALNI_H_V_1(T):
   return -65638
def LLA_AL_ALNI_H_V_2(T):
   return 0
def LLA_AL_ALNI_V_H_0(T):
   return -2945766
def LLA_AL_ALNI_V_H_1(T):
   return -3001368
def LLA_AL_ALNI_V_H_2(T):
   return 0
def LLA_AL_ALNI_V_V_0(T):
   return -190916
def LLA_AL_ALNI_V_V_1(T):
   return 5910
def LLA_AL_ALNI_V_V_2(T):
   return 0
def LLA_NI_ALNI_H_H_0(T):
   return -519565
def LLA_NI_ALNI_H_H_1(T):
   return -836669
def LLA_NI_ALNI_H_H_2(T):
   return 0
def LLA_NI_ALNI_H_V_0(T):
   return 574142
def LLA_NI_ALNI_H_V_1(T):
   return 3457689
def LLA_NI_ALNI_H_V_2(T):
   return 0
def LLA_NI_ALNI_V_V_0(T):
   return -218243
def LLA_NI_ALNI_V_V_1(T):
   return -99998
def LLA_NI_ALNI_V_V_2(T):
   return 0
