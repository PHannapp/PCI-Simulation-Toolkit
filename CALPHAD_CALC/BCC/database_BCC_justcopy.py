import numpy as np
# reference -------------------------------------------------------------------------------
# Al
def GHSERAL(T):
    return -7976.15+137.093038*T-24.3671976*T*np.log(T)-1.884662E-3*T**2-0.877664E-6*T**3+74092*T**(-1)
def GBCCAL(T):
    return 10083-4.813*T+GHSERAL(T)
# Co
def GHSERCO(T):
    return  310.241+133.36601*T-25.0861*T*np.log(T)-2.654739E-3*T**2-0.17348E-6*T**3+72527*T**(-1)
def GBCCCO(T):
    return 2938-0.7138*T+GHSERCO(T)
# Cr
def GHSERCR(T):
    return  -8856.94 + 157.48*T - 26.908*T*np.log(T) + 1.89435E-3*T**2 - 1.47721E-6*T**3 + 139250*T**(-1)
# Fe
def GHSERFE(T):
    return 1225.7+124.134*T-23.5143*T*np.log(T)-4.39752E-3*T**2-0.058927E-6*T**3+77359*T**(-1)
# Hf
def GHSERHF(T):
    return  -6987.297+110.744026*T-22.7075*T*np.log(T)-4.146145E-3*T**2-0.000477E-6*T**3-22590*T**(-1)
def GBCCHF(T):
    return  5370.703+103.836026*T-22.8995*T*np.log(T)-4.206605E-3*T**2+0.871923E-6*T**3-22590*T**(-1)-0.1446E-9*T**4
# Mg
def GHSERMG(T):
    return  -8367.34+143.675547*T-26.1849782*T*np.log(T)+0.4858E-3*T**2-1.393669E-6*T**3+78950*T**(-1)
def GBCCMG(T):
    return 3100-2.1*T+GHSERMG(T)
# Mn
def GBCCMN(T):
    return  -3235.3+127.85*T-23.7*T*np.log(T)-7.44271E-3*T**2+60000*T**(-1)
# Mo
def GHSERMO(T):
    return  -7746.302+131.9197*T-23.56414*T*np.log(T)-3.443396E-3*T**2+0.566283E-6*T**3+65812*T**(-1)-0.130927E-9*T**4
# Nb
def GHSERNB(T):
    return  -8519.353 + 142.045475*T - 26.4711*T*np.log(T) + 0.203475E-3*T**2 - 0.35012E-6*T**3 + 93399*T**(-1)
# Ni
def GHSERNI(T):
    return -5179.159+117.854*T-22.096*T*np.log(T)-0.0048407*T**2
def GBCCNI(T):
    return 8715.084-3.556*T+GHSERNI(T)
# Sc
def GHSERSC(T):
    return -8689.547+153.48097*T-28.1882*T*np.log(T)+3.21892E-3*T**2-1.64531E-6*T**3+72177*T**(-1)
def GBCCSC(T):
    return -6709.819+152.456835*T-28.1882*T*np.log(T)+3.21892E-3*T**2-1.64531E-6*T**3+72177*T**(-1)
# Ta
def GHSERTA(T):
    return -7285.889+119.139857*T-23.7592624*T*np.log(T)-2.623033E-3*T**2+0.170109E-6*T**3-3293*T**(-1)
# Ti
def GBCCTI(T):
    return -1272.064 + 134.71418*T - 25.5768*T*np.log(T) - 0.663845E-3*T**2 - 0.278803E-6*T**3 + 7208*T**(-1)
# V
def GHSERV(T):
    return -7930.43 + 133.346053*T - 24.134*T*np.log(T) - 3.098E-3*T**2 + 0.12175E-6*T**3 + 69460*T**(-1)
# Zr
def GBCCZR(T):
    return -525.539 + 124.9457*T - 25.607406*T*np.log(T) - 0.340084E-3*T**2 - 0.009729E-6*T**3 + 25233*T**(-1) - 0.076143E-9*T**4
# H2
def GHSERH2(T):
    return -9522.97393+78.5273873*T-31.35707*T*np.log(T)+0.0027589925*T**2-0.000000746390667*T**3+56582.3*T**(-1)

# Endmembers -------------------------------------------------------------------------------
# Al-H, Kim, 2012, trimmed
def GAL_V(T):
 return GBCCAL(T) 
def GAL_H(T):
 return GHSERAL(T) + 0.5*GHSERH2(T) + 14491.72646478753 + 95.72990611652753*T
# Co-H, Kim, 2012, trimmed
def GCO_V(T):
 return GBCCCO(T) 
def GCO_H(T):
 return GHSERCO(T) + 0.5*GHSERH2(T) + 43380.76551470249 + 99.14590617365882*T
# Cr-H, Joubert, 2022, trimmed
def GCR_V(T):
 return GHSERCR(T) 
def GCR_H(T):
 return GHSERCR(T) + 0.5*GHSERH2(T) + 39388 + 40.99*T
# Fe-H, Zinkevich, 2002
def GFE_V(T):
 return GHSERFE(T)
def GFE_H(T):
 return GHSERFE(T) + 0.5*GHSERH2(T) + 171098.3995 - 986.034062*T + 1155.942278*T*np.log(T) - 0.0894166963*T**2
# Hf-H, from Ti
def GHF_V(T):
 return GBCCHF(T) 
def GHF_H(T):
 return GHSERHF(T) + 0.5*GHSERH2(T) - 63853.977074357914 + 53.931669040437335*T
# Mg-H, from Nb
def GMG_V(T):
 return GBCCMG(T) 
def GMG_H(T):
 return GHSERMG(T) + 0.5*GHSERH2(T) -37171 + 49.5*T
# Mn-H, Liang, 2018, trimmed
def GMN_V(T):
 return GBCCMN(T)
def GMN_H(T):
 return GBCCMN(T) + 0.5*GHSERH2(T) + 17942.66000005651 + 32.809047881388054*T
# Mo-H, from Cr
def GMO_V(T):
 return GHSERMO(T) 
def GMO_H(T):
 return GHSERMO(T) + 0.5*GHSERH2(T) + 39388 + 40.99*T
# Nb-H, Dottor, 2023, trimmed
def GNB_V(T):
 return GHSERNB(T)
def GNB_H(T):
 return GHSERNB(T) + 0.5*GHSERH2(T) -37171 + 49.5*T
# Ni-H, Bourgeois, 2015, DFT?, from Mn
def GNI_V(T):
 return GBCCNI(T)
def GNI_H(T):
 return GHSERNI(T) + 0.5*GHSERH2(T) + 17942.66000005651 + 32.809047881388054*T
# Sc-H, from Zr - 10,000
def GSC_V(T):
 return GBCCSC(T) 
def GSC_H(T):
 return GHSERSC(T) + 0.5*GHSERH2(T) - 80754 + 49*T
# Ta-H, from Nb
def GTA_V(T):
 return GHSERTA(T)
def GTA_H(T):
 return GHSERTA(T) + 0.5*GHSERH2(T) -37171 + 49.5*T
# Ti-H, Wang, 2010, trimmed
def GTI_V(T):
 return GBCCTI(T)
def GTI_H(T):
 return GBCCTI(T) + 0.5*GHSERH2(T) - 63853.977074357914 + 53.931669040437335*T
# V-H, Ukita, 2008, trimmed
def GVA_V(T):
 return GHSERV(T) 
def GVA_H(T):
 return GHSERV(T) + 0.5*GHSERH2(T) -31868.74827689512 + 57.75093422388118*T
# Zr-H, Dottor, 2022, trimmed
def GZR_V(T):
 return GBCCZR(T)
def GZR_H(T):
 return GBCCZR(T) + 0.5*GHSERH2(T) - 70754 + 49*T

# Interaction parameters -------------------------------------------------------------------
# Al-H, Kim, 2012, trimmed
def LAL_HV_0(T):
 return -32952.15621861388
def LAL_HV_1(T):
 return 3386.324467644827
def LAL_HV_2(T):
 return 2142.9665707247473
# Co-H, Kim, 2012, trimmed
def LCO_HV_0(T):
 return -27395.408491324408
def LCO_HV_1(T):
 return 3387.5993206968546
def LCO_HV_2(T):
 return 2143.627146792412
# Cr-H, Joubert, 2022, trimmed
def LCR_HV_0(T):
 return -249.67 + 4.924*T
# Fe-H, Zinkevich, 2002
# Hf-H, from Ti
def LHF_HV_0(T):
 return 4596.183547575047 + 4.9047784187241605*T
def LHF_HV_1(T):
 return -2199.6983316091896 + 4.3666463928399395*T
# Mg-H, from Nb
def LMG_HV_0(T):
 return 4288.015361005298
def LMG_HV_1(T):
 return -3958.619672869934
def LMG_HV_2(T):
 return 954.7056583190248
def LMG_HV_3(T):
 return 2200.688918342813
# Mn-H, Liang, 2018, trimmed
def LMN_HV_0(T):
 return 6801.616434023487
# Mo-H, from Cr
def LMO_HV_0(T):
 return -249.67 + 4.924*T
# Nb-H, Dottor, 2023, trimmed
def LNB_HV_0(T):
 return 4288.015361005298
def LNB_HV_1(T):
 return -3958.619672869934
def LNB_HV_2(T):
 return 954.7056583190248
def LNB_HV_3(T):
 return 2200.688918342813
# Ni-H, Bourgeois, 2015, DFT?, from Mn
def LNI_HV_0(T):
 return 6801.616434023487
# Sc-H, from , from Zr
def LSC_HV_0(T):
 return 6924.78 + 5.621*T
def LSC_HV_1(T):
 return -9430.8 + 11.9*T
# Ta-H, from Nb
def LTA_HV_0(T):
 return 4288.015361005298
def LTA_HV_1(T):
 return -3958.619672869934
def LTA_HV_2(T):
 return 954.7056583190248
def LTA_HV_3(T):
 return 2200.688918342813
# Ti-H, Wang, 2010, trimmed
def LTI_HV_0(T):
 return 4596.183547575047 + 4.9047784187241605*T
def LTI_HV_1(T):
 return -2199.6983316091896 + 4.3666463928399395*T
# V-H, Ukita, 2008, trimmed
def LVA_HV_0(T):
 return  -3079.67392323499
def LVA_HV_1(T):
 return -2908.9816522488286
def LVA_HV_2(T):
 return 2575.840406187378
# Zr-H, Dottor, 2022, trimmed
def LZR_HV_0(T):
 return 6924.78 + 5.621*T
def LZR_HV_1(T):
 return -9430.8 + 11.9*T


# Metal - Metal - in Metal ------
# Al-Co, Liu, 2016
def LALCO_V_0(T):
  return -138500+34.62*T + 56531-37.04*T
# Al-Cr, Liu, 2016
def LALCR_V_0(T):
  return -54900 + 10*T
# Al-Fe, Ostrowska, 2019
def LALFE_V_0(T):
  return -122960 + 32*T
def LALFE_V_1(T):
  return 2945.2
# Al-Mg, Hallstedt, 2023
def LALMG_V_0(T):
  return +1593+2.149*T
def LALMG_V_1(T):
  return +1014-0.660*T
def LALMG_V_2(T):
  return -673
# Al-Mn, Hallstedt, 2023
def LALMN_V_0(T):
  return -111336.2 + 44.48059*T
def LALMN_V_1(T):
  return -68691.8 + 44.52653*T
# Al-Nb, He,2015
def LALNB_V_0(T):
  return -47470
def LALNB_V_1(T):
  return +77910-42*T
# Al-Ni, Liu, 2016
def LALNi_V_0(T):
  return -152397.3+26.40575*T - 52440.88+11.30117*T
# Al-Ti, Witusiewicz, 2007
def LALTI_V_0(T):
  return -132903 + 39.961*T
def LALTI_V_1(T):
  return 4890
def LALTI_V_2(T):
  return 400
# Al-V, Kroupa, 2017
def LALVA_V_0(T):
  return -89400+16.46*T
def LALVA_V_1(T):
  return -6000
# Co-Cr, Liu, 2016
def LCOCR_V_0(T):
  return +17208-13.519*T
def LCOCR_V_1(T):
  return -5470
# Co-Fe, Ostrowska, 2019
def LCOFE_V_0(T):
  return -20205 + 14.8*T+0.9845*T*np.log(T)-7.6434E-03*T**2
def LCOFE_V_2(T):
  return 1316
# Co-Mn, Hallstedt, 2023
def LCOMN_V_0(T):
  return -23945
# Co-Nb, He, 2015
def LCONB_V_0(T):
  return 7400
# Co-Ni, Liu, 2016
def LCONI_V_0(T):
  return 2000
# Co-Ti, Zhou, 2018
def LCOTI_V_0(T):
  return -92966+12.38*T
# Co-V, Choi, 2019
def LCOVA_V_0(T):
  return -35366.216+17.435*T
# Cr-Fe, Hallstedt, 2023
def LCRFE_V_0(T):
  return +20500-9.68*T
# Cr-Mg, Cui, 2017
def LCRMG_V_0(T):
  return 80*T
# Cr-Mn, Hallstedt, 2023
def LCRMN_V_0(T):
  return -20328+18.7339*T
def LCRMN_V_1(T):
  return -9162+4.4183*T
# Cr-Nb, Peng, 2016
def LCRNB_V_0(T):
    return 62044.08-23.1246*T
def LCRNB_V_1(T):
    return 38777.36-19.31*T
# Cr-Ni, Liu, 2016
def LCRNI_V_0(T):
  return +17170 - 11.8199*T
def LCRNI_V_1(T):
  return +34418 - 11.8577*T
# Cr-Sc, Wang, 2014
def LCRSC_V_0(T):
  return 47900 + 3*T
def LCRSC_V_1(T):
  return 25800
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
# Fe-Mg, Dilner, 2016 -> 
def LFEMG_V_0(T):
 return 65700
#Fe-Mn, Kim, 2015
def LFEMN_V_0(T):
 return -2759 + 1.24*T
# Fe-Nb, Lu, 2017
def LFENB_V_0(T):
 return 1238+5.08*T
def LFENB_V_1(T):
 return 11532-11.33*T
# Fe-Ni, Hallstedt, 2023
def LFENI_V_0(T):
  return -7500
def LFENI_V_1(T):
  return 8500-5*T
# Fe-Ti,  Guo, 2012
def LFETI_V_0(T):
 return -66435+22.512*T
def LFETI_V_1(T):
 return +5968-5.083*T
def LFETI_V_2(T):
 return +31454-20.139*T
#Fe-V, Guo, 2012
def LFEVA_V_0(T):
 return -21427+6.846*T
def LFEVA_V_1(T):
 return +7345-1.509*T
# Fe-Zr, Lu, 2017
def LFEZR_V_0(T):
 return -12200+4.5*T
def LFEZR_V_1(T):
 return -3410
def LFEZR_V_2(T):
 return -15000
# Hf-Mo, Xiao, 2024
def LHFMO_V_0(T):
 return -24872+8.905*T
def LHFMO_V_1(T):
 return 65515-29.900*T
def LHFMO_V_2(T):
 return 2159-5.847*T
# Hf-Nb, Ghosh, 2002
def LHFNB_V_0(T):
 return -5835.23 + 16.32278*T
def LHFNB_V_1(T):
 return -3660.25 + 6.44118*T
# Hf-Ti, Xiao, 2024
def LHFTI_V_0(T):
 return 3003.24-7.41140*T
# Hf-Zr, Bittermann, 2002
def LHFZR_V_0(T):
 return -8386.21
# Mg-Mn, Wang, 2015
def LMGMN_V_0(T):
 return 85000
# Mg-Nb, not available!!!!
# Mg-Ni, Wang, 2015
def LMGNI_V_0(T):
 return 80*T
# Mg-Ti, not available!!!!
# Mg-V, not available!!!!
# Mn-Nb, Khvan, 2012
def LMNNB_V_0(T):
 return 6305.5
# Mn-Ni, Hallstedt, 2023
def LMNNI_V_0(T):
  return -3508.42803-23.7885576*T
# Mn-Ti, Khan, 2016
def LMNTI_V_0(T):
 return -66236.1 + 28.2691*T
def LMNTI_V_1(T):
 return 6000
# Mn-V, Huang, 1991
def LMNVA_V_0(T):
 return -10000 
# Mn-Zr, Flandorfer, 2021
def LMNZR_V_0(T):
 return -54491.27+30.8108*T
def LMNZR_V_1(T):
 return -6418.27+0.51346*T
# Mo-Nb
def LMONB_V_0(T):
 return -68202.6+29.85596*T
def LMONB_V_1(T):
 return +8201.3
# Mo-Ti, Xiao, 2024
def LMOTI_V_0(T):
 return 20249-13.474*T
def LMOTI_V_1(T):
 return -22963+8.598*T
def LMOTI_V_2(T):
 return -8171
# Mo-V, Yang, 2021
def LMOVA_V_0(T):
 return -14800.1+13.0712*T
def LMOVA_V_1(T):
 return 9
def LMOVA_V_2(T):
 return -33315+12.9*T
# Mo-Zr, Che, 2024
def LMOZR_V_0(T):
 return -2768 +13.411*T
def LMOZR_V_1(T):
 return -9152 +6.856*T
# Nb-Ni, Huang, 2019, TC?
def LNBNI_V_0(T):
    return -18724.3 + 5.02405*T
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
# Ni-Sc, Zhu, 2015
def LNISC_V_0(T):
 return -9.6609240E+04 - 1.8291157E+01*T
def LNISC_V_1(T):
 return -1.0866503E+04 - 1.0636766E+01*T
# Ni-Ti, Zhou, 2018
def LNITI_V_0(T):
  return -97427.4+12.112*T
def LNITI_V_1(T):
  return -32315
# Ni-V, Huang, 2016
def LNIVA_V_0(T):
 return -7835.86451 - 5.88123433*T
def LNIVA_V_1(T):
 return +38301.8097-28.0663731*T
# Sc-Ti
# Sc-V
# Ti-V, Cui, 2016
def LTIVA_V_0(T):
 return 6523.17
def LTIVA_V_1(T):
 return 2025.39
# Ti-Zr, Sridar, 2017
def LTIZR_V_0(T):
    return -4100 + 3.47*T
# V-Zr, Cui, 2016
def LVAZR_V_0(T):
 return 17872.99 + 8.7539*T
def LVAZR_V_1(T):
 return -3208.60 + 4.8481*T