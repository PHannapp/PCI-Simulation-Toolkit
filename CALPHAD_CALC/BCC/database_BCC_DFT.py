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
def GHSERFE(T):
    return 1225.7+124.134*T-23.5143*T*np.log(T)-4.39752E-3*T**2-0.058927E-6*T**3+77359*T**(-1)
# Mn
def GBCCMN(T):
    return  -3235.3+127.85*T-23.7*T*np.log(T)-7.44271E-3*T**2+60000*T**(-1)
# Ni
def GHSERNI(T):
    return -5179.159+117.854*T-22.096*T*np.log(T)-0.0048407*T**2
# Ta
def GHSERTA(T):
    return -7285.889+119.139857*T-23.7592624*T*np.log(T)-2.623033E-3*T**2+0.170109E-6*T**3-3293*T**(-1)
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
 return GHSERNB(T)
def GNB_H(T):
 return GHSERNB(T) + 0.5*GHSERH2(T) -37171 + 49.5*T
# Cr-H, Joubert, 2022
def GCR_V(T):
 return GHSERCR(T) 
def GCR_H(T):
 return GHSERCR(T) + 0.5*GHSERH2(T) + 39388 + 40.99*T
# Fe-H, Zinkevich, 2002
def GFE_V(T):
 return GHSERFE(T)
def GFE_H(T):
 return GHSERFE(T) + 0.5*GHSERH2(T) + 171098.3995 - 986.034062*T + 1155.942278*T*np.log(T) - 0.0894166963*T**2
#Mn-H, Liang, 2018
def GMN_V(T):
 return GBCCMN(T)
def GMN_H(T):
 return GBCCMN(T) + 0.5*GHSERH2(T) + 17942.66000005651 + 32.809047881388054*T
# Ni-H, Bourgeois, 2015
def GNI_V(T):
 return GHSERNI(T)
def GNI_H(T):
 return GHSERNI(T) + 0.5*GHSERH2(T) -9.15446340e3 + 6.97721747e1*T - 5.82437764e-3*T**2
# Ta unknown
def GTA_V(T):
 return GHSERTA(T)
def GTA_H(T):
 return GHSERTA(T) + 0.5*GHSERH2(T) - 17624 + 78.88*T
# Ti-H, Wang, 2010
def GTI_V(T):
 return GBCCTI(T)
def GTI_H(T):
 return GBCCTI(T) + 0.5*GHSERH2(T) - 63853.977074357914 + 53.931669040437335*T
# V-H, Ukita, 2008
def GVA_V(T):
 return GHSERV(T) 
def GVA_H(T):
 return GHSERV(T) + 0.5*GHSERH2(T) + 9898.79389341804 - 47.17945618543901*T
# Zr-H, Dottor, 2022
def GZR_V(T):
 return GBCCZR(T)
def GZR_H(T):
 return GBCCZR(T) + 0.5*GHSERH2(T) - 70754 + 49*T

# Interaction parameters -------------------------------------------------------------------
# Nb-H, Dottor, 2023
def LNB_HV_0(T):
 return 4288.015361005298
def LNB_HV_1(T):
 return -3958.619672869934
def LNB_HV_2(T):
 return 954.7056583190248
def LNB_HV_3(T):
 return 2200.688918342813
# Cr-H, Joubert, 2022
def LCR_HV_0(T):
 return -249.67 + 4.924*T
# Fe-H, Zinkevich, 2002
# Mn-H, fake
def LMN_HV_0(T):
 return 6801.616434023487
# Ni-H, Bourgeois, 2015
def LNI_HV_0(T):
 return +1.44616815e4 - 4.45244666*T
def LNI_HV_1(T):
 return -7.00967408e2
# Ti-H, Wang, 2010
def LTI_HV_0(T):
 return 4596.183547575047 + 4.9047784187241605*T
def LTI_HV_1(T):
 return -2199.6983316091896 + 4.3666463928399395*T
# V-H, Ukita, 2008
def LVA_HV_0(T):
 return  -2623.725430155838 - 1.1873375540181603*T
def LVA_HV_1(T):
 return -299.51381017669814 -9.628595925549977*T
def LVA_HV_2(T):
 return -136.34784782571725 + 5.362060640340153*T
# Zr-H, Dottor, 2022
def LZR_HV_0(T):
 return 6924.78 + 5.621*T
def LZR_HV_1(T):
 return -9430.8 + 11.9*T


# Nb-Cr, Lu, 2015
def LNBCR_V_0(T):
    return +43600 - 13.6*T
def LNBCR_V_1(T):
    return +13755.2 - 6.8*T
# Nb-Fe
# Nb-Mn
# Nb-Ni, Huang, 2019, TC?
def LNBNI_V_0(T):
    return -18724.3 + 5.02405*T
# Nb-Ta
# Nb-Ti, Zhang, 2022
def LNBTI_V_0(T):
    return 14000
def LNBTI_V_2(T):
    return 2500
# Nb-V, Kumar, 1994
def LNBVA_V_0(T):
    return +9080
# Nb-Zr, Dottor, 2023
def LNBZR_V_0(T):
    return +14389.575 + 4.2621477*T
def LNBZR_V_1(T):
    return +3417.0688
# Cr-Fe
# Cr-Mn
# Cr-Ni
# Cr-Ta
# Cr-Ti, Gosh, 2001
def LCRTI_V_0(T):
    return -2247.87 + 9.14144*T
def LCRTI_V_1(T):
    return 198.73
# Cr-V, Gosh, 2001
def LCRVA_V_0(T):
    return -8253.85 - 3.61592*T
def LCRVA_V_1(T):
    return 7494.82 - 8.69424*T
def LCRVA_V_2(T):
    return -17599.07 + 10.13142*T
# Cr-Zr, Lu, 2015
def LCRZR_V_0(T):
    return +43600 - 13.6*T
def LCRZR_V_1(T):
    return +13755.2 - 6.8*T
# Fe-Mn, Kim, 2015
def LFEMN_V_0(T):
 return -2759 + 1.24*T
# Fe-Ni
# Fe-Ta
# Fe-Ti,  Guo, 2012
def LFETI_V_0(T):
 return -66435+22.512*T
def LFETI_V_1(T):
 return +5968-5.083*T
def LFETI_V_2(T):
 return +31454-20.139*T
# Fe-V, Guo, 2012
def LFEVA_V_0(T):
 return -21427+6.846*T
def LFEVA_V_1(T):
 return +7345-1.509*T
# Fe-Zr
# Mn-Ni
# Mn-Ta
# Mn-Ti, Khan, 2016
def LMNTI_V_0(T):
 return -66236.1 + 28.2691*T
def LMNTI_V_1(T):
 return 6000
# Mn-V, Huang, 1991
def LMNVA_V_0(T):
 return -10000
# Mn-Zr
# Ni-Ta
# Ni-Ti
# Ni-V, Huang, 2016
def LNIVA_V_0(T):
 return -7835.86451 - 5.88123433*T
def LNIVA_V_1(T):
 return +38301.8097-28.0663731*T
# Ni-Zr
# Ti-V, Cui, 2016
def LTIVA_V_0(T):
 return 6523.17
def LTIVA_V_1(T):
 return 2025.39
# Ti-Zr, Sridar, 2017
def LTIZR_V_0(T):
    return -4100 + 3.47*T
# Ta-Ti
# Ta-V
# Ta-Zr
# Ti-V
# Ti-Zr
# V-Zr, Cui, 2016
def LVAZR_V_0(T):
 return 17872.99 + 8.7539*T
def LVAZR_V_1(T):
 return -3208.60 + 4.8481*T

