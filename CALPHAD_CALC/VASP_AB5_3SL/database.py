import numpy as np
# reference -------------------------------------------------------------------------------
def GHSERFE(T):
    return 1225.7 + 124.134*T - 23.5143*T*np.log(T) - 0.00439752*T**2 - 5.8927E-08*T**3 + 77359*T**(-1)
def GHSERNI(T):
    return -5179.159+117.854*T-22.096*T*np.log(T)-0.0048407*T**2
def GHSERAL(T):
    return -7976.15 + 137.093038*T - 24.3671976*T*np.log(T) - 1.884662E-3*T**2 - 0.877664E-6*T**3 + 74092*T**(-1)
def GHSERH2(T):
    return -9522.97393+78.5273873*T-31.35707*T*np.log(T)+0.0027589925*T**2-0.000000746390667*T**3+56582.3*T**(-1)
def GHSERCU(T):
    return -7770.458 + 130.485235 * T - 24.112392 * T * np.log(T) - 2.65684E-3 * T**2 + 0.129223E-6 * T**3 + 52478 * T**(-1)
def GHSERLA(T):
    return -7968.403 + 120.284604*T - 26.34*T*np.log(T) - 0.001295165*T**2
def GHSERCA(T):
    return -4955.062 + 72.794266 * T - 16.3138 * T * np.log(T) - 11.10455E-3 * T**2 - 133574 * T**(-1)
def GHSERCE(T):
    return -7160.519+84.23022*T-22.3664*T*np.log(T)-0.0067103*T**2-0.000000320773*T**3-18117*T**(-1)

# Endmembers -------------------------------------------------------------------------------
def GCE_AL_H_H_H(T):
   return 1*GHSERCE(T) + 5*GHSERAL(T) + 3.5*GHSERH2(T) + (178534) + (427.0 * T)
def GCE_AL_H_H_V(T):
   return 1*GHSERCE(T) + 5*GHSERAL(T) + 2.0*GHSERH2(T) + (37881) + (244.0 * T)
def GCE_AL_H_V_H(T):
   return 1*GHSERCE(T) + 5*GHSERAL(T) + 2.0*GHSERH2(T) + (19942) + (244.0 * T)
def GCE_AL_H_V_V(T):
   return 1*GHSERCE(T) + 5*GHSERAL(T) + 0.5*GHSERH2(T) + (-45267) + (61.0 * T)
def GCE_AL_V_H_V(T):
   return 1*GHSERCE(T) + 5*GHSERAL(T) + 1.5*GHSERH2(T) + (68274) + (183.0 * T)
def GCE_AL_V_V_H(T):
   return 1*GHSERCE(T) + 5*GHSERAL(T) + 1.5*GHSERH2(T) + (-31669) + (183.0 * T)
def GCE_AL_V_V_V(T):
   return 1*GHSERCE(T) + 5*GHSERAL(T) + 0*GHSERH2(T) + (-114000) + (0 * T)
def GCE_NI_H_H_H(T):
   return 1*GHSERCE(T) + 5*GHSERNI(T) + 3.5*GHSERH2(T) + (-298685) + (427.0 * T)
def GCE_NI_H_H_V(T):
   return 1*GHSERCE(T) + 5*GHSERNI(T) + 2.0*GHSERH2(T) + (-231278) + (244.0 * T)
def GCE_NI_H_V_H(T):
   return 1*GHSERCE(T) + 5*GHSERNI(T) + 2.0*GHSERH2(T) + (-248480) + (244.0 * T)
def GCE_NI_H_V_V(T):
   return 1*GHSERCE(T) + 5*GHSERNI(T) + 0.5*GHSERH2(T) + (-219878) + (61.0 * T)
def GCE_NI_V_H_H(T):
   return 1*GHSERCE(T) + 5*GHSERNI(T) + 3.0*GHSERH2(T) + (-173598) + (366.0 * T)
def GCE_FE_H_H_H(T):
   return 1*GHSERCE(T) + 5*GHSERFE(T) + 3.5*GHSERH2(T) + (-3969665) + (427.0 * T)
def GCE_FE_H_H_V(T):
   return 1*GHSERCE(T) + 5*GHSERFE(T) + 2.0*GHSERH2(T) + (-3948592) + (244.0 * T)
def GCE_FE_H_V_H(T):
   return 1*GHSERCE(T) + 5*GHSERFE(T) + 2.0*GHSERH2(T) + (-3954689) + (244.0 * T)
def GCE_FE_H_V_V(T):
   return 1*GHSERCE(T) + 5*GHSERFE(T) + 0.5*GHSERH2(T) + (-3892603) + (61.0 * T)
def GCE_FE_V_H_H(T):
   return 1*GHSERCE(T) + 5*GHSERFE(T) + 3.0*GHSERH2(T) + (-3839098) + (366.0 * T)
def GCE_FE_V_H_V(T):
   return 1*GHSERCE(T) + 5*GHSERFE(T) + 1.5*GHSERH2(T) + (-3785243) + (183.0 * T)
def GCE_FE_V_V_H(T):
   return 1*GHSERCE(T) + 5*GHSERFE(T) + 1.5*GHSERH2(T) + (-3938396) + (183.0 * T)
def GCE_FE_V_V_V(T):
   return 1*GHSERCE(T) + 5*GHSERFE(T) + 0*GHSERH2(T) + (-3817240) + (0 * T)
def GCE_CU_H_H_H(T):
   return 1*GHSERCE(T) + 5*GHSERCU(T) + 3.5*GHSERH2(T) + (-32416) + (427.0 * T)
def GCE_CU_H_V_H(T):
   return 1*GHSERCE(T) + 5*GHSERCU(T) + 2.0*GHSERH2(T) + (-97402) + (244.0 * T)
def GCE_CU_H_V_V(T):
   return 1*GHSERCE(T) + 5*GHSERCU(T) + 0.5*GHSERH2(T) + (-41270) + (61.0 * T)
def GLA_AL_H_H_H(T):
   return 1*GHSERLA(T) + 5*GHSERAL(T) + 3.5*GHSERH2(T) + (154557) + (427.0 * T)
def GLA_AL_H_H_V(T):
   return 1*GHSERLA(T) + 5*GHSERAL(T) + 2.0*GHSERH2(T) + (54561) + (244.0 * T)
def GLA_AL_H_V_H(T):
   return 1*GHSERLA(T) + 5*GHSERAL(T) + 2.0*GHSERH2(T) + (-11035) + (244.0 * T)
def GLA_AL_H_V_V(T):
   return 1*GHSERLA(T) + 5*GHSERAL(T) + 0.5*GHSERH2(T) + (-72959) + (61.0 * T)
def GLA_AL_V_H_H(T):
   return 1*GHSERLA(T) + 5*GHSERAL(T) + 3.0*GHSERH2(T) + (108260) + (366.0 * T)
def GLA_AL_V_H_V(T):
   return 1*GHSERLA(T) + 5*GHSERAL(T) + 1.5*GHSERH2(T) + (33957) + (183.0 * T)
def GLA_AL_V_V_H(T):
   return 1*GHSERLA(T) + 5*GHSERAL(T) + 1.5*GHSERH2(T) + (-57781) + (183.0 * T)
def GLA_AL_V_V_V(T):
   return 1*GHSERLA(T) + 5*GHSERAL(T) + 0*GHSERH2(T) + (-130942) + (0 * T)
def GLA_NI_H_H_H(T):
   return 1*GHSERLA(T) + 5*GHSERNI(T) + 3.5*GHSERH2(T) + (-318867) + (427.0 * T)
def GLA_NI_H_H_V(T):
   return 1*GHSERLA(T) + 5*GHSERNI(T) + 2.0*GHSERH2(T) + (-247596) + (244.0 * T)
def GLA_NI_H_V_H(T):
   return 1*GHSERLA(T) + 5*GHSERNI(T) + 2.0*GHSERH2(T) + (-231307) + (244.0 * T)
def GLA_NI_H_V_V(T):
   return 1*GHSERLA(T) + 5*GHSERNI(T) + 0.5*GHSERH2(T) + (-189129) + (61.0 * T)
def GLA_NI_V_H_H(T):
   return 1*GHSERLA(T) + 5*GHSERNI(T) + 3.0*GHSERH2(T) + (-177388) + (366.0 * T)
def GLA_NI_V_H_V(T):
   return 1*GHSERLA(T) + 5*GHSERNI(T) + 1.5*GHSERH2(T) + (-221057) + (183.0 * T)
def GLA_NI_V_V_H(T):
   return 1*GHSERLA(T) + 5*GHSERNI(T) + 1.5*GHSERH2(T) + (-219938) + (183.0 * T)
def GLA_NI_V_V_V(T):
   return 1*GHSERLA(T) + 5*GHSERNI(T) + 0*GHSERH2(T) + (-91638) + (0 * T)
def GLA_FE_H_H_H(T):
   return 1*GHSERLA(T) + 5*GHSERFE(T) + 3.5*GHSERH2(T) + (-3949190) + (427.0 * T)
def GLA_FE_H_H_V(T):
   return 1*GHSERLA(T) + 5*GHSERFE(T) + 2.0*GHSERH2(T) + (-3905929) + (244.0 * T)
def GLA_FE_H_V_H(T):
   return 1*GHSERLA(T) + 5*GHSERFE(T) + 2.0*GHSERH2(T) + (-3883361) + (244.0 * T)
def GLA_FE_H_V_V(T):
   return 1*GHSERLA(T) + 5*GHSERFE(T) + 0.5*GHSERH2(T) + (-3811893) + (61.0 * T)
def GLA_FE_V_H_H(T):
   return 1*GHSERLA(T) + 5*GHSERFE(T) + 3.0*GHSERH2(T) + (-3863481) + (366.0 * T)
def GLA_FE_V_H_V(T):
   return 1*GHSERLA(T) + 5*GHSERFE(T) + 1.5*GHSERH2(T) + (-3770673) + (183.0 * T)
def GLA_FE_V_V_H(T):
   return 1*GHSERLA(T) + 5*GHSERFE(T) + 1.5*GHSERH2(T) + (-3879171) + (183.0 * T)
def GLA_FE_V_V_V(T):
   return 1*GHSERLA(T) + 5*GHSERFE(T) + 0*GHSERH2(T) + (-3739520) + (0 * T)
def GLA_CU_H_H_H(T):
   return 1*GHSERLA(T) + 5*GHSERCU(T) + 3.5*GHSERH2(T) + (-85805) + (427.0 * T)
def GLA_CU_H_H_V(T):
   return 1*GHSERLA(T) + 5*GHSERCU(T) + 2.0*GHSERH2(T) + (-60104) + (244.0 * T)
def GLA_CU_H_V_H(T):
   return 1*GHSERLA(T) + 5*GHSERCU(T) + 2.0*GHSERH2(T) + (-120312) + (244.0 * T)
def GLA_CU_H_V_V(T):
   return 1*GHSERLA(T) + 5*GHSERCU(T) + 0.5*GHSERH2(T) + (-72804) + (61.0 * T)
def GLA_CU_V_H_H(T):
   return 1*GHSERLA(T) + 5*GHSERCU(T) + 3.0*GHSERH2(T) + (-67164) + (366.0 * T)
def GLA_CU_V_H_V(T):
   return 1*GHSERLA(T) + 5*GHSERCU(T) + 1.5*GHSERH2(T) + (-25756) + (183.0 * T)
def GLA_CU_V_V_H(T):
   return 1*GHSERLA(T) + 5*GHSERCU(T) + 1.5*GHSERH2(T) + (-128419) + (183.0 * T)
def GLA_CU_V_V_V(T):
   return 1*GHSERLA(T) + 5*GHSERCU(T) + 0*GHSERH2(T) + (-73094) + (0 * T)
def GCA_AL_H_H_H(T):
   return 1*GHSERCA(T) + 5*GHSERAL(T) + 3.5*GHSERH2(T) + (219845) + (427.0 * T)
def GCA_AL_H_H_V(T):
   return 1*GHSERCA(T) + 5*GHSERAL(T) + 2.0*GHSERH2(T) + (122784) + (244.0 * T)
def GCA_AL_H_V_H(T):
   return 1*GHSERCA(T) + 5*GHSERAL(T) + 2.0*GHSERH2(T) + (59110) + (244.0 * T)
def GCA_AL_H_V_V(T):
   return 1*GHSERCA(T) + 5*GHSERAL(T) + 0.5*GHSERH2(T) + (9914) + (61.0 * T)
def GCA_AL_V_H_H(T):
   return 1*GHSERCA(T) + 5*GHSERAL(T) + 3.0*GHSERH2(T) + (178822) + (366.0 * T)
def GCA_AL_V_H_V(T):
   return 1*GHSERCA(T) + 5*GHSERAL(T) + 1.5*GHSERH2(T) + (37253) + (183.0 * T)
def GCA_AL_V_V_H(T):
   return 1*GHSERCA(T) + 5*GHSERAL(T) + 1.5*GHSERH2(T) + (57083) + (183.0 * T)
def GCA_AL_V_V_V(T):
   return 1*GHSERCA(T) + 5*GHSERAL(T) + 0*GHSERH2(T) + (-62889) + (0 * T)
def GCA_NI_H_H_H(T):
   return 1*GHSERCA(T) + 5*GHSERNI(T) + 3.5*GHSERH2(T) + (-241928) + (427.0 * T)
def GCA_NI_H_H_V(T):
   return 1*GHSERCA(T) + 5*GHSERNI(T) + 2.0*GHSERH2(T) + (-174264) + (244.0 * T)
def GCA_NI_H_V_H(T):
   return 1*GHSERCA(T) + 5*GHSERNI(T) + 2.0*GHSERH2(T) + (-147900) + (244.0 * T)
def GCA_NI_H_V_V(T):
   return 1*GHSERCA(T) + 5*GHSERNI(T) + 0.5*GHSERH2(T) + (-84504) + (61.0 * T)
def GCA_NI_V_H_H(T):
   return 1*GHSERCA(T) + 5*GHSERNI(T) + 3.0*GHSERH2(T) + (-103811) + (366.0 * T)
def GCA_NI_V_H_V(T):
   return 1*GHSERCA(T) + 5*GHSERNI(T) + 1.5*GHSERH2(T) + (-113947) + (183.0 * T)
def GCA_NI_V_V_H(T):
   return 1*GHSERCA(T) + 5*GHSERNI(T) + 1.5*GHSERH2(T) + (-134479) + (183.0 * T)
def GCA_NI_V_V_V(T):
   return 1*GHSERCA(T) + 5*GHSERNI(T) + 0*GHSERH2(T) + (14227) + (0 * T)
def GCA_FE_H_H_H(T):
   return 1*GHSERCA(T) + 5*GHSERFE(T) + 3.5*GHSERH2(T) + (-3861068) + (427.0 * T)
def GCA_FE_H_H_V(T):
   return 1*GHSERCA(T) + 5*GHSERFE(T) + 2.0*GHSERH2(T) + (-3825482) + (244.0 * T)
def GCA_FE_H_V_H(T):
   return 1*GHSERCA(T) + 5*GHSERFE(T) + 2.0*GHSERH2(T) + (-3803882) + (244.0 * T)
def GCA_FE_H_V_V(T):
   return 1*GHSERCA(T) + 5*GHSERFE(T) + 0.5*GHSERH2(T) + (-3759084) + (61.0 * T)
def GCA_FE_V_H_H(T):
   return 1*GHSERCA(T) + 5*GHSERFE(T) + 3.0*GHSERH2(T) + (-3753799) + (366.0 * T)
def GCA_FE_V_H_V(T):
   return 1*GHSERCA(T) + 5*GHSERFE(T) + 1.5*GHSERH2(T) + (-3667185) + (183.0 * T)
def GCA_FE_V_V_H(T):
   return 1*GHSERCA(T) + 5*GHSERFE(T) + 1.5*GHSERH2(T) + (-3823363) + (183.0 * T)
def GCA_FE_V_V_V(T):
   return 1*GHSERCA(T) + 5*GHSERFE(T) + 0*GHSERH2(T) + (-3697693) + (0 * T)
def GCA_CU_H_H_H(T):
   return 1*GHSERCA(T) + 5*GHSERCU(T) + 3.5*GHSERH2(T) + (-59781) + (427.0 * T)
def GCA_CU_H_H_V(T):
   return 1*GHSERCA(T) + 5*GHSERCU(T) + 2.0*GHSERH2(T) + (-27759) + (244.0 * T)
def GCA_CU_H_V_H(T):
   return 1*GHSERCA(T) + 5*GHSERCU(T) + 2.0*GHSERH2(T) + (-34604) + (244.0 * T)
def GCA_CU_H_V_V(T):
   return 1*GHSERCA(T) + 5*GHSERCU(T) + 0.5*GHSERH2(T) + (-39681) + (61.0 * T)
def GCA_CU_V_H_H(T):
   return 1*GHSERCA(T) + 5*GHSERCU(T) + 3.0*GHSERH2(T) + (31141) + (366.0 * T)
def GCA_CU_V_H_V(T):
   return 1*GHSERCA(T) + 5*GHSERCU(T) + 1.5*GHSERH2(T) + (23676) + (183.0 * T)
def GCA_CU_V_V_H(T):
   return 1*GHSERCA(T) + 5*GHSERCU(T) + 1.5*GHSERH2(T) + (-58801) + (183.0 * T)
def GCA_CU_V_V_V(T):
   return 1*GHSERCA(T) + 5*GHSERCU(T) + 0*GHSERH2(T) + (-40947) + (0 * T)

# Interaction parameters -------------------------------------------------------------------
