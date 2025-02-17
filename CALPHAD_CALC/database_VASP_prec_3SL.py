import numpy as np
# reference -------------------------------------------------------------------------------
def GHSERCE(T):
    return -7160.519+84.23022*T-22.3664*T*np.log(T)-0.0067103*T**2-0.000000320773*T**3-18117*T**(-1)
def GHSERAL(T):
    return -7976.15 + 137.093038*T - 24.3671976*T*np.log(T) - 1.884662E-3*T**2 - 0.877664E-6*T**3 + 74092*T**(-1)
def GHSERNI(T):
    return -5179.159+117.854*T-22.096*T*np.log(T)-0.0048407*T**2
def GHSERLA(T):
    return -7968.403 + 120.284604*T - 26.34*T*np.log(T) - 0.001295165*T**2
def GHSERH2(T):
    return -9522.97393+78.5273873*T-31.35707*T*np.log(T)+0.0027589925*T**2-0.000000746390667*T**3+56582.3*T**(-1)

# Endmembers -------------------------------------------------------------------------------
def GCE_AL_H_H_H(T):
   return 1*GHSERCE(T) + 5*GHSERAL(T) + 3.5*GHSERH2(T) + (170237) + (427.0 * T)
def GCE_AL_H_H_V(T):
   return 1*GHSERCE(T) + 5*GHSERAL(T) + 2.0*GHSERH2(T) + (165729) + (244.0 * T)
def GCE_AL_H_V_H(T):
   return 1*GHSERCE(T) + 5*GHSERAL(T) + 2.0*GHSERH2(T) + (60724) + (244.0 * T)
def GCE_AL_H_V_V(T):
   return 1*GHSERCE(T) + 5*GHSERAL(T) + 0.5*GHSERH2(T) + (2322) + (61.0 * T)
def GCE_AL_V_H_H(T):
   return 1*GHSERCE(T) + 5*GHSERAL(T) + 3.0*GHSERH2(T) + (143667) + (366.0 * T)
def GCE_AL_V_H_V(T):
   return 1*GHSERCE(T) + 5*GHSERAL(T) + 1.5*GHSERH2(T) + (116763) + (183.0 * T)
def GCE_AL_V_V_H(T):
   return 1*GHSERCE(T) + 5*GHSERAL(T) + 1.5*GHSERH2(T) + (-53057) + (183.0 * T)
def GCE_AL_V_V_V(T):
   return 1*GHSERCE(T) + 5*GHSERAL(T) + 0*GHSERH2(T) + (-129060) + (0 * T)
def GCE_NI_H_H_H(T):
   return 1*GHSERCE(T) + 5*GHSERNI(T) + 3.5*GHSERH2(T) + (-271898) + (427.0 * T)
def GCE_NI_H_H_V(T):
   return 1*GHSERCE(T) + 5*GHSERNI(T) + 2.0*GHSERH2(T) + (-195423) + (244.0 * T)
def GCE_NI_H_V_H(T):
   return 1*GHSERCE(T) + 5*GHSERNI(T) + 2.0*GHSERH2(T) + (-210412) + (244.0 * T)
def GCE_NI_H_V_V(T):
   return 1*GHSERCE(T) + 5*GHSERNI(T) + 0.5*GHSERH2(T) + (-184045) + (61.0 * T)
def GCE_NI_V_H_H(T):
   return 1*GHSERCE(T) + 5*GHSERNI(T) + 3.0*GHSERH2(T) + (-248338) + (366.0 * T)
def GCE_NI_V_H_V(T):
   return 1*GHSERCE(T) + 5*GHSERNI(T) + 1.5*GHSERH2(T) + (-183163) + (183.0 * T)
def GCE_NI_V_V_H(T):
   return 1*GHSERCE(T) + 5*GHSERNI(T) + 1.5*GHSERH2(T) + (-217652) + (183.0 * T)
def GCE_NI_V_V_V(T):
   return 1*GHSERCE(T) + 5*GHSERNI(T) + 0*GHSERH2(T) + (-201441) + (0 * T)
def GLA_AL_H_H_H(T):
   return 1*GHSERLA(T) + 5*GHSERAL(T) + 3.5*GHSERH2(T) + (146626) + (427.0 * T)
def GLA_AL_H_H_V(T):
   return 1*GHSERLA(T) + 5*GHSERAL(T) + 2.0*GHSERH2(T) + (119370) + (244.0 * T)
def GLA_AL_H_V_H(T):
   return 1*GHSERLA(T) + 5*GHSERAL(T) + 2.0*GHSERH2(T) + (31331) + (244.0 * T)
def GLA_AL_H_V_V(T):
   return 1*GHSERLA(T) + 5*GHSERAL(T) + 0.5*GHSERH2(T) + (-29877) + (61.0 * T)
def GLA_AL_V_H_H(T):
   return 1*GHSERLA(T) + 5*GHSERAL(T) + 3.0*GHSERH2(T) + (117973) + (366.0 * T)
def GLA_AL_V_H_V(T):
   return 1*GHSERLA(T) + 5*GHSERAL(T) + 1.5*GHSERH2(T) + (78206) + (183.0 * T)
def GLA_AL_V_V_H(T):
   return 1*GHSERLA(T) + 5*GHSERAL(T) + 1.5*GHSERH2(T) + (-76800) + (183.0 * T)
def GLA_AL_V_V_V(T):
   return 1*GHSERLA(T) + 5*GHSERAL(T) + 0*GHSERH2(T) + (-151821) + (0 * T)
def GLA_NI_H_H_H(T):
   return 1*GHSERLA(T) + 5*GHSERNI(T) + 3.5*GHSERH2(T) + (-287977) + (427.0 * T)
def GLA_NI_H_H_V(T):
   return 1*GHSERLA(T) + 5*GHSERNI(T) + 2.0*GHSERH2(T) + (-220253) + (244.0 * T)
def GLA_NI_H_V_H(T):
   return 1*GHSERLA(T) + 5*GHSERNI(T) + 2.0*GHSERH2(T) + (-196931) + (244.0 * T)
def GLA_NI_H_V_V(T):
   return 1*GHSERLA(T) + 5*GHSERNI(T) + 0.5*GHSERH2(T) + (-167639) + (61.0 * T)
def GLA_NI_V_H_H(T):
   return 1*GHSERLA(T) + 5*GHSERNI(T) + 3.0*GHSERH2(T) + (-259141) + (366.0 * T)
def GLA_NI_V_H_V(T):
   return 1*GHSERLA(T) + 5*GHSERNI(T) + 1.5*GHSERH2(T) + (-189095) + (183.0 * T)
def GLA_NI_V_V_H(T):
   return 1*GHSERLA(T) + 5*GHSERNI(T) + 1.5*GHSERH2(T) + (-185285) + (183.0 * T)
def GLA_NI_V_V_V(T):
   return 1*GHSERLA(T) + 5*GHSERNI(T) + 0*GHSERH2(T) + (-165893) + (0 * T)

# Interaction parameters -------------------------------------------------------------------
