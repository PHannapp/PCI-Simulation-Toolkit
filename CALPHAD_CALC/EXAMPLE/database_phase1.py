import numpy as np
# reference -------------------------------------------------------------------------------
# Nb
def GHSERNB(T):
    return  -8519.353 + 142.045475*T - 26.4711*T*np.log(T) + 0.203475E-3*T**2 - 0.35012E-6*T**3 + 93399*T**(-1)
def GFCCNB(T):
    return  13500 + 1.7*T + GHSERNB(T)
# Zr
def GHSERZR(T):
    return -7827.595 + 125.64905*T - 24.1618*T*np.log(T) - 4.37791E-3*T**2 + 34971*T**(-1)
def GBCCZR(T):
    return -525.539 + 124.9457*T - 25.607406*T*np.log(T) - 0.340084E-3*T**2 - 0.009729E-6*T**3 + 25233*T**(-1) - 0.076143E-9*T**4
def GFCCZR(T):
    return 7600 - 0.9*T + GHSERZR(T)
# H2
def GHSERH2(T):
    return -9522.97393+78.5273873*T-31.35707*T*np.log(T)+0.0027589925*T**2-0.000000746390667*T**3+56582.3*T**(-1)

# Endmembers -------------------------------------------------------------------------------
# Nb-H, Dottor, 2023
def GNB_V(T):
 return GHSERNB(T)
def GNB_H(T):
 return GHSERNB(T) + 0.5*GHSERH2(T) -37171 + 49.5*T
# Zr-H, Dottor, 2022
def GZR_V(T):
 return GBCCZR(T)
def GZR_H(T):
 return GBCCZR(T) + 0.5*GHSERH2(T) - 70754 + 49*T
# Interaction parameters -------------------------------------------------------------------
# Nb-H, Dottor, 2023
def LNB_HV_0(T):
 return 5990.47
def LNB_HV_1(T):
 return -3112.46
def LNB_HV_2(T):
 return 1848.64
def LNB_HV_3(T):
 return 2928.25
# Zr-H, Dottor, 2022
def LZR_HV_0(T):
 return 6924.78 + 5.621*T
def LZR_HV_1(T):
 return -9430.8 + 11.9*T
# Nb-Zr, Dottor, 2023
def LNBZR_V_0(T):
    return +14389.575 + 4.2621477*T
def LNBZR_V_1(T):
    return +3417.0688