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
   return 1*GHSERCE(T) + 5*GHSERAL(T) + 3.5*GHSERH2(T) + (178913) + (427.0 * T)
def GCE_AL_H_H_V(T):
   return 1*GHSERCE(T) + 5*GHSERAL(T) + 2.0*GHSERH2(T) + (169216) + (244.0 * T)
def GCE_AL_H_V_H(T):
   return 1*GHSERCE(T) + 5*GHSERAL(T) + 2.0*GHSERH2(T) + (71900) + (244.0 * T)
def GCE_AL_H_V_V(T):
   return 1*GHSERCE(T) + 5*GHSERAL(T) + 0.5*GHSERH2(T) + (11533) + (61.0 * T)
def GCE_AL_V_H_H(T):
   return 1*GHSERCE(T) + 5*GHSERAL(T) + 3.0*GHSERH2(T) + (154461) + (366.0 * T)
def GCE_AL_V_H_V(T):
   return 1*GHSERCE(T) + 5*GHSERAL(T) + 1.5*GHSERH2(T) + (118610) + (183.0 * T)
def GCE_AL_V_V_H(T):
   return 1*GHSERCE(T) + 5*GHSERAL(T) + 1.5*GHSERH2(T) + (-40742) + (183.0 * T)
def GCE_AL_V_V_V(T):
   return 1*GHSERCE(T) + 5*GHSERAL(T) + 0*GHSERH2(T) + (-119121) + (0 * T)
def GCE_NI_H_H_H(T):
   return 1*GHSERCE(T) + 5*GHSERNI(T) + 3.5*GHSERH2(T) + (-274438) + (427.0 * T)
def GCE_NI_H_H_V(T):
   return 1*GHSERCE(T) + 5*GHSERNI(T) + 2.0*GHSERH2(T) + (-194384) + (244.0 * T)
def GCE_NI_H_V_H(T):
   return 1*GHSERCE(T) + 5*GHSERNI(T) + 2.0*GHSERH2(T) + (-213977) + (244.0 * T)
def GCE_NI_H_V_V(T):
   return 1*GHSERCE(T) + 5*GHSERNI(T) + 0.5*GHSERH2(T) + (-185709) + (61.0 * T)
def GCE_NI_V_H_H(T):
   return 1*GHSERCE(T) + 5*GHSERNI(T) + 3.0*GHSERH2(T) + (-251901) + (366.0 * T)
def GCE_NI_V_H_V(T):
   return 1*GHSERCE(T) + 5*GHSERNI(T) + 1.5*GHSERH2(T) + (-185923) + (183.0 * T)
def GCE_NI_V_V_H(T):
   return 1*GHSERCE(T) + 5*GHSERNI(T) + 1.5*GHSERH2(T) + (-219239) + (183.0 * T)
def GCE_NI_V_V_V(T):
   return 1*GHSERCE(T) + 5*GHSERNI(T) + 0*GHSERH2(T) + (-202280) + (0 * T)
def GLA_AL_H_H_H(T):
   return 1*GHSERLA(T) + 5*GHSERAL(T) + 3.5*GHSERH2(T) + (155113) + (427.0 * T)
def GLA_AL_H_H_V(T):
   return 1*GHSERLA(T) + 5*GHSERAL(T) + 2.0*GHSERH2(T) + (123107) + (244.0 * T)
def GLA_AL_H_V_H(T):
   return 1*GHSERLA(T) + 5*GHSERAL(T) + 2.0*GHSERH2(T) + (46309) + (244.0 * T)
def GLA_AL_H_V_V(T):
   return 1*GHSERLA(T) + 5*GHSERAL(T) + 0.5*GHSERH2(T) + (-19328) + (61.0 * T)
def GLA_AL_V_H_H(T):
   return 1*GHSERLA(T) + 5*GHSERAL(T) + 3.0*GHSERH2(T) + (131944) + (366.0 * T)
def GLA_AL_V_H_V(T):
   return 1*GHSERLA(T) + 5*GHSERAL(T) + 1.5*GHSERH2(T) + (80824) + (183.0 * T)
def GLA_AL_V_V_H(T):
   return 1*GHSERLA(T) + 5*GHSERAL(T) + 1.5*GHSERH2(T) + (-61331) + (183.0 * T)
def GLA_AL_V_V_V(T):
   return 1*GHSERLA(T) + 5*GHSERAL(T) + 0*GHSERH2(T) + (-146205) + (0 * T)
def GLA_NI_H_H_H(T):
   return 1*GHSERLA(T) + 5*GHSERNI(T) + 3.5*GHSERH2(T) + (-294242) + (427.0 * T)
def GLA_NI_H_H_V(T):
   return 1*GHSERLA(T) + 5*GHSERNI(T) + 2.0*GHSERH2(T) + (-223252) + (244.0 * T)
def GLA_NI_H_V_H(T):
   return 1*GHSERLA(T) + 5*GHSERNI(T) + 2.0*GHSERH2(T) + (-199098) + (244.0 * T)
def GLA_NI_H_V_V(T):
   return 1*GHSERLA(T) + 5*GHSERNI(T) + 0.5*GHSERH2(T) + (-169786) + (61.0 * T)
def GLA_NI_V_H_H(T):
   return 1*GHSERLA(T) + 5*GHSERNI(T) + 3.0*GHSERH2(T) + (-260587) + (366.0 * T)
def GLA_NI_V_H_V(T):
   return 1*GHSERLA(T) + 5*GHSERNI(T) + 1.5*GHSERH2(T) + (-193559) + (183.0 * T)
def GLA_NI_V_V_H(T):
   return 1*GHSERLA(T) + 5*GHSERNI(T) + 1.5*GHSERH2(T) + (-185352) + (183.0 * T)
def GLA_NI_V_V_V(T):
   return 1*GHSERLA(T) + 5*GHSERNI(T) + 0*GHSERH2(T) + (-165723) + (0 * T)

# Interaction parameters -------------------------------------------------------------------
