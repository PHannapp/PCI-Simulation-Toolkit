import numpy as np
# reference -------------------------------------------------------------------------------
def GHSERLA(T):
    return -7968.403 + 120.284604*T - 26.34*T*np.log(T) - 0.001295165*T**2
def GHSERAL(T):
    return -7976.15 + 137.093038*T - 24.3671976*T*np.log(T) - 1.884662E-3*T**2 - 0.877664E-6*T**3 + 74092*T**(-1)
def GHSERNI(T):
    return -5179.159+117.854*T-22.096*T*np.log(T)-0.0048407*T**2
def GHSERCE(T):
    return -7160.519+84.23022*T-22.3664*T*np.log(T)-0.0067103*T**2-0.000000320773*T**3-18117*T**(-1)
def GHSERH2(T):
    return -9522.97393+78.5273873*T-31.35707*T*np.log(T)+0.0027589925*T**2-0.000000746390667*T**3+56582.3*T**(-1)

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
   return -53174.92307692308
def LLA_ALNI_H_H_1(T):
   return -19155.923076923078
def LLA_ALNI_H_H_2(T):
   return -62188.769230769234
def LLA_ALNI_H_V_0(T):
   return -114050.0
def LLA_ALNI_H_V_1(T):
   return -47593.57142857143
def LLA_ALNI_H_V_2(T):
   return -50088.57142857143
def LLA_ALNI_V_H_0(T):
   return -70507.16666666667
def LLA_ALNI_V_H_1(T):
   return -7406.333333333333
def LLA_ALNI_V_H_2(T):
   return -22234.75
def LLA_AL_HV_H_0(T):
   return -27852.384615384617
def LLA_AL_HV_H_1(T):
   return -4426.384615384615
def LLA_AL_HV_H_2(T):
   return -32127.923076923078
def LLA_AL_HV_V_0(T):
   return 3605.714285714286
def LLA_AL_HV_V_1(T):
   return -668.2857142857143
def LLA_AL_HV_V_2(T):
   return 276.85714285714283


def LLA_AL_H_HV_0(T):
   return -26807.46153846154
def LLA_AL_H_HV_1(T):
   return -2795.3076923076924
def LLA_AL_H_HV_2(T):
   return -73443.0
def LLA_AL_V_HV_0(T):
   return 6226.333333333333
def LLA_AL_V_HV_1(T):
   return -51606.833333333336
def LLA_AL_V_HV_2(T):
   return -35398.416666666664


def LLA_NI_HV_H_0(T):
   return 2247
def LLA_NI_HV_H_1(T):
   return -1079
def LLA_NI_HV_H_2(T):
   return 2832
def LLA_NI_HV_V_0(T):
   return 10772
def LLA_NI_HV_V_1(T):
   return 2540
def LLA_NI_HV_V_2(T):
   return 6870

def LLA_NI_H_HV_0(T):
   return 108701
def LLA_NI_H_HV_1(T):
   return -23481
def LLA_NI_H_HV_2(T):
   return 22427
def LLA_NI_V_HV_0(T):
   return 107048
def LLA_NI_V_HV_1(T):
   return -104275
def LLA_NI_V_HV_2(T):
   return 67849