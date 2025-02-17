import numpy as np
# reference -------------------------------------------------------------------------------
def GHSERLA(T):
    return -7968.403 + 120.284604*T - 26.34*T*np.log(T) - 0.001295165*T**2
def GHSERNI(T):
    return -5179.159+117.854*T-22.096*T*np.log(T)-0.0048407*T**2
def GHSERAL(T):
    return -7976.15 + 137.093038*T - 24.3671976*T*np.log(T) - 1.884662E-3*T**2 - 0.877664E-6*T**3 + 74092*T**(-1)
def GHSERFE(T):
    return 1225.7 + 124.134*T - 23.5143*T*np.log(T) - 0.00439752*T**2 - 5.8927E-08*T**3 + 77359*T**(-1)
def GHSERCU(T):
    return -7770.458 + 130.485235 * T - 24.112392 * T * np.log(T) - 2.65684E-3 * T**2 + 0.129223E-6 * T**3 + 52478 * T**(-1)
def GHSERCA(T):
    return -4955.062 + 72.794266 * T - 16.3138 * T * np.log(T) - 11.10455E-3 * T**2 - 133574 * T**(-1)
def GHSERH2(T):
    return -9522.97393+78.5273873*T-31.35707*T*np.log(T)+0.0027589925*T**2-0.000000746390667*T**3+56582.3*T**(-1)
def GHSERCE(T):
    return -7160.519+84.23022*T-22.3664*T*np.log(T)-0.0067103*T**2-0.000000320773*T**3-18117*T**(-1)

# Endmembers -------------------------------------------------------------------------------
def GCE_AL_H_H(T):
   return 1*GHSERCE(T) + 5*GHSERAL(T) + 3.5*GHSERH2(T) + (179060) + (427.0 * T)
def GCE_AL_H_V(T):
   return 1*GHSERCE(T) + 5*GHSERAL(T) + 0.5*GHSERH2(T) + (-45496) + (61.0 * T)
def GCE_NI_H_H(T):
   return 1*GHSERCE(T) + 5*GHSERNI(T) + 3.5*GHSERH2(T) + (-298685) + (427.0 * T)
def GCE_FE_H_H(T):
   return 1*GHSERCE(T) + 5*GHSERFE(T) + 3.5*GHSERH2(T) + (-3969678) + (427.0 * T)
def GCE_FE_H_V(T):
   return 1*GHSERCE(T) + 5*GHSERFE(T) + 0.5*GHSERH2(T) + (-3892666) + (61.0 * T)
def GCE_FE_V_H(T):
   return 1*GHSERCE(T) + 5*GHSERFE(T) + 3.0*GHSERH2(T) + (-3839167) + (366.0 * T)
def GCE_CU_H_H(T):
   return 1*GHSERCE(T) + 5*GHSERCU(T) + 3.5*GHSERH2(T) + (-32416) + (427.0 * T)
def GCE_CU_H_V(T):
   return 1*GHSERCE(T) + 5*GHSERCU(T) + 0.5*GHSERH2(T) + (-41527) + (61.0 * T)
def GLA_AL_H_H(T):
   return 1*GHSERLA(T) + 5*GHSERAL(T) + 3.5*GHSERH2(T) + (154558) + (427.0 * T)
def GLA_AL_H_V(T):
   return 1*GHSERLA(T) + 5*GHSERAL(T) + 0.5*GHSERH2(T) + (-73093) + (61.0 * T)
def GLA_AL_V_H(T):
   return 1*GHSERLA(T) + 5*GHSERAL(T) + 3.0*GHSERH2(T) + (108505) + (366.0 * T)
def GLA_AL_V_V(T):
   return 1*GHSERLA(T) + 5*GHSERAL(T) + 0*GHSERH2(T) + (-130875) + (0 * T)

def GLA_NI_H_H(T):
   return 1*GHSERLA(T) + 5*GHSERNI(T) + 3.5*GHSERH2(T) + (-318867) + (427.0 * T)
def GLA_NI_H_V(T):
   return 1*GHSERLA(T) + 5*GHSERNI(T) + 0.5*GHSERH2(T) + (-189153) + (61.0 * T)
def GLA_NI_V_H(T):
   return 1*GHSERLA(T) + 5*GHSERNI(T) + 3.0*GHSERH2(T) + (-177694) + (366.0 * T)
def GLA_NI_V_V(T):
   return 1*GHSERLA(T) + 5*GHSERNI(T) + 0*GHSERH2(T) + (-89522) + (0 * T)

def GLA_FE_H_H(T):
   return 1*GHSERLA(T) + 5*GHSERFE(T) + 3.5*GHSERH2(T) + (-3949113) + (427.0 * T)
def GLA_FE_H_V(T):
   return 1*GHSERLA(T) + 5*GHSERFE(T) + 0.5*GHSERH2(T) + (-3811943) + (61.0 * T)
def GLA_FE_V_H(T):
   return 1*GHSERLA(T) + 5*GHSERFE(T) + 3.0*GHSERH2(T) + (-3863412) + (366.0 * T)
def GLA_FE_V_V(T):
   return 1*GHSERLA(T) + 5*GHSERFE(T) + 0*GHSERH2(T) + (-3739481) + (0 * T)
def GLA_CU_H_H(T):
   return 1*GHSERLA(T) + 5*GHSERCU(T) + 3.5*GHSERH2(T) + (-85805) + (427.0 * T)
def GLA_CU_H_V(T):
   return 1*GHSERLA(T) + 5*GHSERCU(T) + 0.5*GHSERH2(T) + (-73376) + (61.0 * T)
def GLA_CU_V_H(T):
   return 1*GHSERLA(T) + 5*GHSERCU(T) + 3.0*GHSERH2(T) + (-67135) + (366.0 * T)
def GLA_CU_V_V(T):
   return 1*GHSERLA(T) + 5*GHSERCU(T) + 0*GHSERH2(T) + (-72993) + (0 * T)
def GCA_AL_H_H(T):
   return 1*GHSERCA(T) + 5*GHSERAL(T) + 3.5*GHSERH2(T) + (219845) + (427.0 * T)
def GCA_AL_H_V(T):
   return 1*GHSERCA(T) + 5*GHSERAL(T) + 0.5*GHSERH2(T) + (9914) + (61.0 * T)
def GCA_AL_V_H(T):
   return 1*GHSERCA(T) + 5*GHSERAL(T) + 3.0*GHSERH2(T) + (178822) + (366.0 * T)
def GCA_AL_V_V(T):
   return 1*GHSERCA(T) + 5*GHSERAL(T) + 0*GHSERH2(T) + (-62889) + (0 * T)
def GCA_NI_H_H(T):
   return 1*GHSERCA(T) + 5*GHSERNI(T) + 3.5*GHSERH2(T) + (-241928) + (427.0 * T)
def GCA_NI_H_V(T):
   return 1*GHSERCA(T) + 5*GHSERNI(T) + 0.5*GHSERH2(T) + (-84536) + (61.0 * T)
def GCA_NI_V_H(T):
   return 1*GHSERCA(T) + 5*GHSERNI(T) + 3.0*GHSERH2(T) + (-103851) + (366.0 * T)
def GCA_NI_V_V(T):
   return 1*GHSERCA(T) + 5*GHSERNI(T) + 0*GHSERH2(T) + (14228) + (0 * T)
def GCA_FE_H_H(T):
   return 1*GHSERCA(T) + 5*GHSERFE(T) + 3.5*GHSERH2(T) + (-3861062) + (427.0 * T)
def GCA_FE_H_V(T):
   return 1*GHSERCA(T) + 5*GHSERFE(T) + 0.5*GHSERH2(T) + (-3759085) + (61.0 * T)
def GCA_FE_V_H(T):
   return 1*GHSERCA(T) + 5*GHSERFE(T) + 3.0*GHSERH2(T) + (-3753800) + (366.0 * T)
def GCA_FE_V_V(T):
   return 1*GHSERCA(T) + 5*GHSERFE(T) + 0*GHSERH2(T) + (-3697691) + (0 * T)
def GCA_CU_H_H(T):
   return 1*GHSERCA(T) + 5*GHSERCU(T) + 3.5*GHSERH2(T) + (-59781) + (427.0 * T)
def GCA_CU_H_V(T):
   return 1*GHSERCA(T) + 5*GHSERCU(T) + 0.5*GHSERH2(T) + (-39682) + (61.0 * T)
def GCA_CU_V_H(T):
   return 1*GHSERCA(T) + 5*GHSERCU(T) + 3.0*GHSERH2(T) + (31141) + (366.0 * T)
def GCA_CU_V_V(T):
   return 1*GHSERCA(T) + 5*GHSERCU(T) + 0*GHSERH2(T) + (-40910) + (0 * T)

# Interaction parameters -------------------------------------------------------------------
