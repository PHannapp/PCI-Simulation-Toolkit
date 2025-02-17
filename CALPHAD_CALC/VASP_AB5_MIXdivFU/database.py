import numpy as np
# reference -------------------------------------------------------------------------------
def GHSERNI(T):
    return -5179.159+117.854*T-22.096*T*np.log(T)-0.0048407*T**2
def GHSERLA(T):
    return -7968.403 + 120.284604*T - 26.34*T*np.log(T) - 0.001295165*T**2
def GHSERH2(T):
    return -9522.97393+78.5273873*T-31.35707*T*np.log(T)+0.0027589925*T**2-0.000000746390667*T**3+56582.3*T**(-1)
def GHSERAL(T):
    return -7976.15 + 137.093038*T - 24.3671976*T*np.log(T) - 1.884662E-3*T**2 - 0.877664E-6*T**3 + 74092*T**(-1)

# Endmembers -------------------------------------------------------------------------------
def GLA_AL_H_H(T):
   return 1*GHSERLA(T) + 5*GHSERAL(T) + 3.5*GHSERH2(T) + (164916) + (427.0 * T)
def GLA_AL_H_V(T):
   return 1*GHSERLA(T) + 5*GHSERAL(T) + 0.5*GHSERH2(T) + (-68843) + (61.0 * T)
def GLA_AL_V_H(T):
   return 1*GHSERLA(T) + 5*GHSERAL(T) + 3.0*GHSERH2(T) + (187045) + (366.0 * T)
def GLA_AL_V_V(T):
   return 1*GHSERLA(T) + 5*GHSERAL(T) + 0*GHSERH2(T) + (-140522) + (0 * T)
def GLA_NI_H_H(T):
   return 1*GHSERLA(T) + 5*GHSERNI(T) + 3.5*GHSERH2(T) + (-288540) + (427.0 * T)
def GLA_NI_H_V(T):
   return 1*GHSERLA(T) + 5*GHSERNI(T) + 0.5*GHSERH2(T) + (-166615) + (61.0 * T)
def GLA_NI_V_H(T):
   return 1*GHSERLA(T) + 5*GHSERNI(T) + 3.0*GHSERH2(T) + (-259020) + (366.0 * T)
def GLA_NI_V_V(T):
   return 1*GHSERLA(T) + 5*GHSERNI(T) + 0*GHSERH2(T) + (-165532) + (0 * T)

# Interaction parameters -------------------------------------------------------------------
def LLA_ALNI_H_H_0(T):
   return -50908.92307692308
def LLA_ALNI_H_H_1(T):
   return -15823.076923076924
def LLA_ALNI_H_H_2(T):
   return -70340.07692307692
def LLA_ALNI_H_V_0(T):
   return -93157.42857142857
def LLA_ALNI_H_V_1(T):
   return -11760.57142857143
def LLA_ALNI_H_V_2(T):
   return 122676.28571428571
def LLA_ALNI_V_H_0(T):
   return -69175.08333333333
def LLA_ALNI_V_H_1(T):
   return -60434.5
def LLA_ALNI_V_H_2(T):
   return -29566.75
def LLA_ALNI_V_V_0(T):
   return -79708.0
def LLA_ALNI_V_V_1(T):
   return -6264.333333333333
def LLA_ALNI_V_V_2(T):
   return -13271.666666666666
def LLA_AL_HV_H_0(T):
   return -37295.230769230766
def LLA_AL_HV_H_1(T):
   return 12157.0
def LLA_AL_HV_H_2(T):
   return -73585.92307692308
def LLA_AL_HV_V_0(T):
   return 7529.857142857143
def LLA_AL_HV_V_1(T):
   return -41131.57142857143
def LLA_AL_HV_V_2(T):
   return -23049.428571428572
def LLA_NI_HV_H_0(T):
   return -454.6923076923077
def LLA_NI_HV_H_1(T):
   return -792.2307692307693
def LLA_NI_HV_H_2(T):
   return -1461.5384615384614
def LLA_NI_HV_V_0(T):
   return 3015.714285714286
def LLA_NI_HV_V_1(T):
   return -32.42857142857143
def LLA_NI_HV_V_2(T):
   return -2426.1428571428573
def LLA_AL_H_HV_0(T):
   return -15343.307692307691
def LLA_AL_H_HV_1(T):
   return -19590.923076923078
def LLA_AL_H_HV_2(T):
   return -13494.615384615385
def LLA_AL_V_HV_0(T):
   return -3973.25
def LLA_AL_V_HV_1(T):
   return -57813.0
def LLA_AL_V_HV_2(T):
   return -54155.916666666664
def LLA_NI_H_HV_0(T):
   return 6761.538461538462
def LLA_NI_H_HV_1(T):
   return -4377.076923076923
def LLA_NI_H_HV_2(T):
   return -6631.692307692308
def LLA_NI_V_HV_0(T):
   return 11905.916666666666
def LLA_NI_V_HV_1(T):
   return -9687.0
def LLA_NI_V_HV_2(T):
   return -620.4166666666666
