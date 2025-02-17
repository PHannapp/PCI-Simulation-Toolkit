import numpy as np
# reference -------------------------------------------------------------------------------
# Nb
def GHSERNB(T):
    return  -8519.353 + 142.045475*T - 26.4711*T*np.log(T) + 0.203475E-3*T**2 - 0.35012E-6*T**3 + 93399*T**(-1)
def GFCCNB(T):
    return  13500 + 1.7*T + GHSERNB(T)
# Cr
def GHSERCR(T):
    return  -8856.94 + 157.48*T - 26.908*T*np.log(T) + 1.89435E-3*T**2 - 1.47721E-6*T**3 + 139250*T**(-1)
def GFCCCR(T):
    return 7284 + 0.163*T + GHSERCR(T)
# Fe
def GFCCFE(T):
    return -236.7+132.416*T-24.6643*T*np.log(T)-3.75752E-3*T**2-0.058927E-6*T**3+77359*T**(-1)
# Mn
def GFCCMN(T):
    return  -3439.3+131.884*T-24.5177*T*np.log(T)-6E-3*T**2+69600*T**(-1)
# Ni
def GHSERNI(T):
    return -5179.159+117.854*T-22.096*T*np.log(T)-0.0048407*T**2
# Ta
def GHSERTA(T):
    return -7285.889+119.139857*T-23.7592624*T*np.log(T)-2.623033E-3*T**2+0.170109E-6*T**3-3293*T**(-1)
def GFCCTA(T):
    return 16000+1.7*T+GHSERTA(T)

# Ti
def GHSERTI(T):
    return -8059.921 + 133.615208*T - 23.9933*T*np.log(T) - 4.777975E-3*T**2 + 0.106716E-6*T**3 + 72636*T**(-1)
def GBCCTI(T):
    return -1272.064 + 134.71418*T - 25.5768*T*np.log(T) - 0.663845E-3*T**2 - 0.278803E-6*T**3 + 7208*T**(-1)
def GFCCTI(T):
    return 6000 - 0.1*T + GHSERTI(T)
# V
def GHSERV(T):
    return -7930.43 + 133.346053*T - 24.134*T*np.log(T) - 3.098E-3*T**2 + 0.12175E-6*T**3 + 69460*T**(-1)
def GFCCV(T):
    return 7500+1.7*T+GHSERV(T)
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
 return GFCCNB(T)
def GNB_H(T):
 return GHSERNB(T) + 1*GHSERH2(T) -63103.503 + 128.56*T
# Cr-H, Joubert, 2022
def GCR_V(T):
 return GFCCCR(T) 
def GCR_H(T):
 return GHSERCR(T) + 1*GHSERH2(T) -16523.433664038977 + 2*6.0127e1 * T + 190.28011257029962*T
#Fe-H, Zinkevich, 2002
def GFE_V(T):
 return GFCCFE(T)
def GFE_H(T):
 return GFCCFE(T) + 1*GHSERH2(T) + 2*(25746.7209 + 48.1250165*T)
#Mn-H, Liang, 2018
def GMN_V(T):
 return GFCCMN(T)
def GMN_H(T):
 return GFCCMN(T) + 1*GHSERH2(T) +4*(-3814.98 + 20.6553*T)
# Ta unknown
def GTA_V(T):
 return GFCCTA(T)
def GTA_H(T):
 return GFCCTA(T) + 1*GHSERH2(T) + 3*(-18400 + 45.7*T)
# Ti-H, Wang, 2010
def GTI_V(T):
 return GFCCTI(T)
def GTI_H(T):
 return GFCCTI(T) + 1*GHSERH2(T) -141748 + 138.601*T
# V-H, Ukita, 2008
def GVA_V(T):
 return GFCCV(T) 
def GVA_H(T):
 return GHSERV(T) + 1*GHSERH2(T) + 3*(-18400 + 45.7*T)
# Zr-H, Dottor, 2022
def GZR_V(T):
 return GFCCZR(T)
def GZR_H(T):
 return GHSERZR(T) + 1*GHSERH2(T) -170490 + 208.2 * T - 9.47*T*np.log(T)
# Interaction parameters -------------------------------------------------------------------
# Nb-H
# Cr-H, fitted
def LCR_HV_0(T):
 return -11610.724907045284 -335.803691010752*T
def LCR_HV_1(T):
 return 443.9871695819126 -267.9854962564192*T
def LCR_HV_2(T):
 return 11182.934529765545 - 124.24172348893339*T
# Fe-H
#Mn-H, Liang, 2018
def LMN_HV_0(T):
 return 2*(8612.59 + 2.2585*T)
# Ta-H
# Ti-H, Wang, 2010
def LTI_HV_0(T):
 return 8229.73
def LTI_HV_1(T):
 return -17080.7
# V-H, Ukita, 2008
def LVA_HV_0(T):
 return 3*(700)
# Zr-H, Dottor, 2022
def LZR_HV_0(T):
 return +14385 - 6*T
def LZR_HV_1(T):
 return -106445 + 87.3*T

# Nb-Cr
# Nb-Fe
# Nb-Mn
# Nb-Ni, Huang, 2019
def LNBNI_V_0(T):
 return -70007.4 - 7.39665*T
def LNBNI_V_1(T):
 return +96115 - 23.07497*T
# Nb-Ta
# Nb-Ti, Zhang, 2022
def LNBTI_V_0(T):
 return 13600
def LNBTI_V_2(T):
 return 2500
# Nb-V
# Nb-Zr
# Cr-Fe
# Cr-Mn
# Cr-Ni
# Cr-Ta
# Cr-Ti
# Cr-V
# Cr-Zr
# Fe-Mn, Kim, 2015
def LFEMN_V_0(T):
 return -7762 + 3.87*T 
def LFEMN_V_1(T):
 return -259
# Fe-Ni
# Fe-Ta
# Fe-Ti,  Guo, 2012
def LFETI_V_0(T):
 return -56022+8.356*T
def LFETI_V_1(T):
 return +4773-4.029*T
def LFETI_V_2(T):
 return +30021-12.614*T
# Fe-V, Guo, 2012
def LFEVA_V_0(T):
 return -15260+1.765*T
# Fe-Zr
# Mn-Ni
# Mn-Ta
# Mn-Ti, Khan, 2016
def LMNTI_V_0(T):
 return -22500
# Mn-V, Huang, 1991
def LMNVA_V_0(T):
 return -11820
# Mn-Zr
# Ni-Ta
# Ni-Ti
# Ni-V, Huang, 2016
def LNIVA_V_0(T):
 return -42546.1601 + 4.60410753*T
def LNIVA_V_1(T):
 return -66881.7021 + 36.4608175*T
def LNIVA_V_2(T):
 return +36212.3605 - 25.5189889*T
# Ni-Zr
# Ta-Ti
# Ta-V
# Ta-Zr
# Ti-V
# Ti-Zr
# V-Zr