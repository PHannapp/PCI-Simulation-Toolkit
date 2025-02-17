import numpy as np
# reference -------------------------------------------------------------------------------
def GHSERAL(T):
    return -7976.15 + 137.093038*T - 24.3671976*T*np.log(T) - 1.884662E-3*T**2 - 0.877664E-6*T**3 + 74092*T**(-1)
def GHSERCE(T):
    return -7160.519+84.23022*T-22.3664*T*np.log(T)-0.0067103*T**2-0.000000320773*T**3-18117*T**(-1)
def GHSERLA(T):
    return -7968.403 + 120.284604*T - 26.34*T*np.log(T) - 0.001295165*T**2
def GHSERH2(T):
    return -9522.97393+78.5273873*T-31.35707*T*np.log(T)+0.0027589925*T**2-0.000000746390667*T**3+56582.3*T**(-1)
def GHSERNI(T):
    return -5179.159+117.854*T-22.096*T*np.log(T)-0.0048407*T**2

# Endmembers -------------------------------------------------------------------------------
def GCE_AL_H_H(T):
   return 1*GHSERCE(T) + 5*GHSERAL(T) + 3.5*GHSERH2(T) + (177478) + (427.0 * T)
def GCE_AL_H_V(T):
   return 1*GHSERCE(T) + 5*GHSERAL(T) + 0.5*GHSERH2(T) + (11522) + (61.0 * T)
def GCE_AL_V_H(T):
   return 1*GHSERCE(T) + 5*GHSERAL(T) + 3.0*GHSERH2(T) + (154435) + (366.0 * T)
def GCE_AL_V_V(T):
   return 1*GHSERCE(T) + 5*GHSERAL(T) + 0*GHSERH2(T) + (-119151) + (0 * T)
def GCE_NI_H_H(T):
   return 1*GHSERCE(T) + 5*GHSERNI(T) + 3.5*GHSERH2(T) + (-274444) + (427.0 * T)
def GCE_NI_H_V(T):
   return 1*GHSERCE(T) + 5*GHSERNI(T) + 0.5*GHSERH2(T) + (-185709) + (61.0 * T)
def GCE_NI_V_H(T):
   return 1*GHSERCE(T) + 5*GHSERNI(T) + 3.0*GHSERH2(T) + (-251901) + (366.0 * T)
def GCE_NI_V_V(T):
   return 1*GHSERCE(T) + 5*GHSERNI(T) + 0*GHSERH2(T) + (-202362) + (0 * T)
def GLA_AL_H_H(T):
   return 1*GHSERLA(T) + 5*GHSERAL(T) + 3.5*GHSERH2(T) + (155114) + (427.0 * T)
def GLA_AL_H_V(T):
   return 1*GHSERLA(T) + 5*GHSERAL(T) + 0.5*GHSERH2(T) + (-19328) + (61.0 * T)
def GLA_AL_V_H(T):
   return 1*GHSERLA(T) + 5*GHSERAL(T) + 3.0*GHSERH2(T) + (131944) + (366.0 * T)
def GLA_AL_V_V(T):
   return 1*GHSERLA(T) + 5*GHSERAL(T) + 0*GHSERH2(T) + (-146205) + (0 * T)
def GLA_NI_H_H(T):
   return 1*GHSERLA(T) + 5*GHSERNI(T) + 3.5*GHSERH2(T) + (-294242) + (427.0 * T)
def GLA_NI_H_V(T):
   return 1*GHSERLA(T) + 5*GHSERNI(T) + 0.5*GHSERH2(T) + (-169786) + (61.0 * T)
def GLA_NI_V_H(T):
   return 1*GHSERLA(T) + 5*GHSERNI(T) + 3.0*GHSERH2(T) + (-260587) + (366.0 * T)
def GLA_NI_V_V(T):
   return 1*GHSERLA(T) + 5*GHSERNI(T) + 0*GHSERH2(T) + (-165689) + (0 * T)

# Interaction parameters -------------------------------------------------------------------
def LLA_ALNI_H_H_0(T):
   return -694604
def LLA_ALNI_H_H_1(T):
   return -252909
def LLA_ALNI_H_H_2(T):
   return -798581
def LLA_ALNI_H_V_0(T):
   return -796808
def LLA_ALNI_H_V_1(T):
   return -329019
def LLA_ALNI_H_V_2(T):
   return -347724
def LLA_ALNI_V_H_0(T):
   return -825006
def LLA_ALNI_V_H_1(T):
   return -189458
def LLA_ALNI_V_H_2(T):
   return -536741
def LLA_ALNI_V_V_0(T):
   return -464614
def LLA_ALNI_V_V_1(T):
   return 1687
def LLA_ALNI_V_V_2(T):
   return -84195


