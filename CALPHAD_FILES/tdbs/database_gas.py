import numpy as np

def F11(p):
    return -1.8692e-09*p
def F12(p):
    return -2.3753e-10*p
def F13(p):
    return -2.50627e-11*p
def F14(p):
    return -3.4483e-08*p
def F15(p):
    return -1.2469e-08*p
def EXPO1(p):
    return np.exp(F11(p))
def EXPO2(p):
    return np.exp(F12(p))
def EXPO3(p):
    return np.exp(F13(p))
def EXPO4(p):
    return np.exp(F14(p))
def EXPO5(p):
    return np.exp(F15(p))
def FF1(p):
    return 5.35e+08*EXPO1(p)
def FF2(p):
    return 4.21e+09*EXPO2(p)
def FF3(p):
    return 3.99e+10*EXPO3(p)
def FF4(p):
    return 29000000*EXPO4(p)
def FF5(p):
    return 80200000*EXPO5(p)

# reference -------------------------------------------------------------------------------
def GHSERH2(T):
    return -9522.97393+78.5273873*T-31.35707*T*np.log(T)+0.0027589925*T**2-0.000000746390667*T**3+56582.3*T**(-1)