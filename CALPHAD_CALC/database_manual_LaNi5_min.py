import numpy as np
# reference -------------------------------------------------------------------------------
def GHSERCE(T):
    return -7160.519+84.23022*T-22.3664*T*np.log(T)-0.0067103*T**2-0.000000320773*T**3-18117*T**(-1)
def GHSERLA(T):
    return -7968.403 + 120.284604*T - 26.34*T*np.log(T) - 0.001295165*T**2
def GHSERAL(T):
    return -7976.15 + 137.093038*T - 24.3671976*T*np.log(T) - 1.884662E-3*T**2 - 0.877664E-6*T**3 + 74092*T**(-1)
def GHSERH2(T):
    return -9522.97393+78.5273873*T-31.35707*T*np.log(T)+0.0027589925*T**2-0.000000746390667*T**3+56582.3*T**(-1)
def GHSERNI(T):
    return -5179.159+117.854*T-22.096*T*np.log(T)-0.0048407*T**2

# Endmembers -------------------------------------------------------------------------------
def GLA_NI_V(T):
   return 4*GHSERLA(T) + 20*GHSERNI(T) + 0*GHSERH2(T) + (-760.5e3) + (0 * T)
def GLA_NI_H(T):
   return 4*GHSERLA(T) + 20*GHSERNI(T) + 31 / 2 *GHSERH2(T) + (-1314.8e3) + (31 / 2 * 130 * T)

# Interaction parameters -------------------------------------------------------------------
def LLA_NI_HV_0(T):
    return 161153.39334951
def LLA_NI_HV_1(T):
    return  -318023.18103531
def LLA_NI_HV_2(T):
    return -205854.78071867

 
  