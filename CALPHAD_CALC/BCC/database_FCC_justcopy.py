import numpy as np
# reference -------------------------------------------------------------------------------
# Al
def GHSERAL(T):
    return -7976.15+137.093038*T-24.3671976*T*np.log(T)-1.884662E-3*T**2-0.877664E-6*T**3+74092*T**(-1)
# Co
def GHSERCO(T):
    return  310.241+133.36601*T-25.0861*T*np.log(T)-2.654739E-3*T**2-0.17348E-6*T**3+72527*T**(-1)
def GFCCCO(T):
    return 427.591-0.615248*T+GHSERCO(T)
# Cr
def GHSERCR(T):
    return  -8856.94 + 157.48*T - 26.908*T*np.log(T) + 1.89435E-3*T**2 - 1.47721E-6*T**3 + 139250*T**(-1)
def GFCCCR(T):
    return 7284 + 0.163*T + GHSERCR(T)
# Fe
def GFCCFE(T):
    return -236.7+132.416*T-24.6643*T*np.log(T)-3.75752E-3*T**2-0.058927E-6*T**3+77359*T**(-1)
# Hf
def GHSERHF(T):
    return  -6987.297+110.744026*T-22.7075*T*np.log(T)-4.146145E-3*T**2-0.000477E-6*T**3-22590*T**(-1)
def GFCCHF(T):
    return 10000-2.2*T+GHSERHF(T)
# Mg
def GHSERMG(T):
    return  -8367.34+143.675547*T-26.1849782*T*np.log(T)+0.4858E-3*T**2-1.393669E-6*T**3+78950*T**(-1)
def GFCCMG(T):
    return 2600-0.9*T+GHSERMG(T)
# Mn
def GFCCMN(T):
    return  -3439.3+131.884*T-24.5177*T*np.log(T)-6E-3*T**2+69600*T**(-1)
# Mo
def GHSERMO(T):
    return  -7746.302+131.9197*T-23.56414*T*np.log(T)-3.443396E-3*T**2+0.566283E-6*T**3+65812*T**(-1)-0.130927E-9*T**4
def GFCCMO(T):
    return 15200+0.63*T+GHSERMO(T)
# Nb
def GHSERNB(T):
    return  -8519.353 + 142.045475*T - 26.4711*T*np.log(T) + 0.203475E-3*T**2 - 0.35012E-6*T**3 + 93399*T**(-1)
def GFCCNB(T):
    return  13500 + 1.7*T + GHSERNB(T)
# Ni
def GHSERNI(T):
    return -5179.159+117.854*T-22.096*T*np.log(T)-0.0048407*T**2
# Sc
def GHSERSC(T):
    return -8689.547+153.48097*T-28.1882*T*np.log(T)+3.21892E-3*T**2-1.64531E-6*T**3+72177*T**(-1)
def GFCCSC(T):
    return 5000+GHSERSC(T)
# Ta
def GHSERTA(T):
    return -7285.889+119.139857*T-23.7592624*T*np.log(T)-2.623033E-3*T**2+0.170109E-6*T**3-3293*T**(-1)
def GFCCTA(T):
    return 16000+1.7*T+GHSERTA(T)
# Ti
def GHSERTI(T):
    return -8059.921 + 133.615208*T - 23.9933*T*np.log(T) - 4.777975E-3*T**2 + 0.106716E-6*T**3 + 72636*T**(-1)
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
def GFCCZR(T):
    return 7600 - 0.9*T + GHSERZR(T)
# H2
def GHSERH2(T):
    return -9522.97393+78.5273873*T-31.35707*T*np.log(T)+0.0027589925*T**2-0.000000746390667*T**3+56582.3*T**(-1)

# Endmembers -------------------------------------------------------------------------------
# Al-H, Qiu, 2004
def GAL_V(T):
 return GHSERAL(T) 
def GAL_H(T):
 return 100000 + GHSERAL(T) + 0.5*GHSERH2(T)
# Co-H, Kim, 2012
def GCO_V(T):
 return GFCCCO(T) 
def GCO_H(T):
 return GHSERCO(T) + 0.5*GHSERH2(T) + 28000 + 70*T
# Cr-H, Joubert, 2022
def GCR_V(T):
 return GFCCCR(T) 
def GCR_H(T):
 return GHSERCR(T) + 0.5*GHSERH2(T) - 7.6176e3 + 6.0127e1 * T
# Fe-H, Zinkevich, 2002
def GFE_V(T):
 return GFCCFE(T)
def GFE_H(T):
 return GFCCFE(T) + 0.5*GHSERH2(T) +25746.7209 + 48.1250165*T
# Hf-H, DFT?
def GHF_V(T):
 return GHSERHF(T) 
# Mg-H, DFT?
def GMG_V(T):
 return GFCCMG(T) 
# Mn-H, Liang, 2018
def GMN_V(T):
 return GFCCMN(T)
def GMN_H(T):
 return GFCCMN(T) + 0.5*GHSERH2(T) +2*(-3814.98 + 20.6553*T)
def LMN_HV_0(T):
 return 2*(8612.59 + 2.2585*T)
# Mo-H, DFT?
def GMO_V(T):
 return GFCCMO(T) 
# Nb-H, Dottor, 2023
def GNB_V(T):
 return GFCCNB(T)
def GNB_H(T):
 return GHSERNB(T) + 1*GHSERH2(T) -63103.503 + 128.56*T
# Ni-H, Bourgeois, 2015
def GNI_V(T):
 return GHSERNI(T)
def GNI_H(T):
 return GHSERNI(T) + 0.5*GHSERH2(T) -9.15446340e3 + 6.97721747e1*T - 5.82437764e-3*T**2
# Sc-H, DFT?
def GSC_V(T):
 return GFCCSC(T) 
# Ta-H, DFT?
def GTA_V(T):
 return GFCCTA(T) 
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
# Hydrogen - Vacancies ------
# Al-H, Qiu, 2004
def LAL_HV_0(T):
 return -45805 + 56.4302*T
# Co-H, Kim, 2012
def LCO_HV_0(T):
 return -25000 - 6.5*T
# Cr-H, Joubert, 2022, extrapolate, nothing? DFT?
# Fe-H, Zinkevich, 2002, extrapolate, nothing? DFT?
# Hf-H, DFT?
# Mg-H, DFT?
# Mn-H
# Mo-H, DFT?
# Nb-H, Joubert, 2022, nothing? DFT?
# Ni-H, Bourgeois, 2015
def LNI_HV_0(T):
 return +1.44616815e4 - 4.45244666*T
def LNI_HV_1(T):
 return -7.00967408e2
# Sc-H, DFT?
# Ta-H, DFT?
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

# Metal - Metal - in Metal ------
# Al-Co, Liu, 2016
def LALCO_V_0(T):
  return -122840 + 22.925*T
def LALCO_V_2(T):
  return +24568 - 4.585*T
# Al-Cr, Liu, 2016
def LALCR_V_0(T):
  return -45900 + 6*T
# Al-Fe, Ostrowska, 2019
def LALFE_V_0(T):
  return -104700 + 30.65*T
def LALFE_V_1(T):
  return +30000-7*T
def LALFE_V_2(T):
  return +32200-17*T
# Al-Mg, Hallstedt, 2023
def LALMG_V_0(T):
  return +1593+2.149*T
def LALMG_V_1(T):
  return +1014-0.660*T
def LALMG_V_2(T):
  return -673
# Al-Mn, Hallstedt, 2023
def LALMN_V_0(T):
  return -69938+27.17958*T
def LALMN_V_1(T):
  return +8248.9
# Al-Nb, He,2015
def LALNB_V_0(T):
  return -83366
def LALNB_V_1(T):
  return -16000
# Al-Ni, Liu, 2016
def LALNI_V_0(T):
  return -162407.75 + 16.212965*T
def LALNI_V_1(T):
  return +73417.798 - 34.914168*T
def LALNI_V_2(T):
  return +33471.014 - 9.8373558*T
def LALNI_V_3(T):
  return -30758.01 + 10.25267*T
# Al-Ti, Witusiewicz, 2007
def LALTI_V_0(T):
  return -119185 + 40.723*T
# Al-V, Kroupa, 2017
def LALVA_V_0(T):
  return -69947.2+12.33*T
# Co-Cr, Liu, 2016
def LCOCR_V_0(T):
  return +1500-9.592*T
# Co-Fe, Ostrowska, 2019
def LCOFE_V_0(T):
  return -9112 + 3.3*T
def LCOFE_V_2(T):
  return 1667
# Co-Mn, Hallstedt, 2023
def LCOMN_V_0(T):
  return -23756
def LCOMN_V_1(T):
  return -2343
# Co-Nb, He, 2015
def LCONB_V_0(T):
  return -34360
def LCONB_V_1(T):
  return +35000 - 22*T
# Co-Ni, Liu, 2016
def LCONI_V_0(T):
  return -800 + 1.2629*T
# Co-Ti, Zhou, 2018
def LCOTI_V_0(T):
  return -77800-7.4*T
def LCOTI_V_1(T):
  return -1300
# Co-V, Choi, 2019
def LCOVA_V_0(T):
  return -57151.224+21.699*T
def LCOVA_V_1(T):
  return -4581.3012
# Cr-Fe, Hallstedt, 2023
def LCRFE_V_0(T):
  return +10833-7.477*T
def LCRFE_V_1(T):
  return 1410
# Cr-Mg, Cui, 2017
def LCRMG_V_0(T):
  return 80*T
# Cr-Mn, Hallstedt, 2023
def LCRMN_V_0(T):
  return -19088+17.5423*T
# Cr-Nb, Peng, 2016, not available
# Cr-Ni, Liu, 2016
def LCRNI_V_0(T):
  return +8030 - 12.8801*T
def LCRNI_V_1(T):
  return +33080 - 16.0362*T
# Cr-Sc, Wang, 2014, not available
# Cr-Ti, tbd?
# Cr-V, tbd?
# Cr-Zr, tbd?
# Fe-Mg, Dilner, 2016
def LFEMG_V_0(T):
 return 65200
#Fe-Mn, Kim, 2015
def LFEMN_V_0(T):
 return -7762 + 3.87*T 
def LFEMN_V_1(T):
 return -259
# Fe-Nb, Lu, 2017
def LFENB_V_0(T):
 return 8144-10.48*T
# Fe-Ni, Hallstedt, 2023
def LFENI_V_0(T):
  return -15500+2.85*T
def LFENI_V_1(T):
  return +14000-4*T
def LFENI_V_2(T):
  return -3000
#Fe-Ti,  Guo, 2012
def LFETI_V_0(T):
 return -56022+8.356*T
def LFETI_V_1(T):
 return +4773-4.029*T
def LFETI_V_2(T):
 return +30021-12.614*T
#Fe-V, Guo, 2012
def LFEVA_V_0(T):
 return -15260+1.765*T
# Fe-Zr, Lu, 2017
def LFEZR_V_0(T):
 return -22000-3*T
def LFEZR_V_1(T):
 return 10000-6.6*T
# Hf-Mo, Xiao, nothing
# Hf-Nb, Ghosh, 2002
def LHFNB_V_0(T):
 return 3887.36 + 0.02358*T
# Hf-Ti, Xiao, 2024, nothing
# Hf-Zr, Bittermann, 2002, nothing
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
 return 16895.03
# Mn-Ni, Hallstedt, 2023
def LMNNI_V_0(T):
  return -20455.7462-11.7930695*T
def LMNNI_V_1(T):
  return +15581.5125-8.49116947*T
#Mn-Ti, Khan, 2016
def LMNTI_V_0(T):
 return -22500
#Mn-V, Huang, 1991
def LMNVA_V_0(T):
 return -11820
# Mn-Zr, Flandorfer, 2021
def LMNZR_V_0(T):
 return +5059.03-0.40472*T
def LMNZR_V_1(T):
 return -26984.11+2.15873*T
# Mo-Nb, Yang, 2021, nothing
# Mo-Ti, Xiao, 2024
def LMOTI_V_0(T):
 return -25505
# Mo-V, Yang, 2021, nothing
# Mo-Zr, Che, 2024, nothing
# Nb-Ni, Huang, 2019
def LNBNI_V_0(T):
 return -70007.4 - 7.39665*T
def LNBNI_V_1(T):
 return +96115 - 23.07497*T
# Nb-Ti, Zhang, 2022
def LNBTI_V_0(T):
 return 13600
def LNBTI_V_2(T):
 return 2500
# Ni-V, Huang, 2016
def LNIVA_V_0(T):
 return -42546.1601 + 4.60410753*T
def LNIVA_V_1(T):
 return -66881.7021 + 36.4608175*T
def LNIVA_V_2(T):
 return +36212.3605 - 25.5189889*T
# Ni-Sc, Zhu, 2015
def LNISC_V_0(T):
 return -1.1308328E+05 - 7.4479178*T
# Ni-Ti, Zhou, 2018
def LNITI_V_0(T):
  return -98143+6.706*T
def LNITI_V_1(T):
  return -62430
# Ni-V, tbd?
# Nb-Zr, tbd?
# Sc-Ti
# Sc-V
# Ti-V, tbd?
# Ti-Zr, tbd?
# V-Zr, tbd?

