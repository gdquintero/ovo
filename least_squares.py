import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def func_zika(t,a,b,c):
    ebt = np.exp(-1.0 * b * t)
    res = (a / b) * t * ebt + (1.0 / b) * ((a / b) - c) * (ebt - 1.0) - c * t
    res = 1.0 - np.exp(res)
    return res

def func_wave(t,a,b,c,d):
    return a * np.exp(b * np.sin(c * t)) + d

def func_original(t,a,b,c,d):
    return a + b * t + c * (t**2) + d * (t**3)

outlier = True
model = 2

if model == 0:
    if outlier:
        df_file = "output/zika_outliers.txt"
    else:
        df_file = "output/zika.txt"

    df = pd.read_csv(df_file,header=None, sep=" ")
    popt, pcov = curve_fit(func_zika,df[0].values,df[1].values,bounds=(0.,np.inf * np.ones(3)))

    with open("output/xstarls.txt",'w') as f:
        f.write("%11.8f\n" % popt[0])
        f.write("%11.8f\n" % popt[1])
        f.write("%11.8f" % popt[2])

elif model == 1:
    df_file = "output/wave.txt"
    df = pd.read_csv(df_file,header=None, sep=" ")
    popt, pcov = curve_fit(func_wave,df[0].values,df[1].values,method='dogbox',maxfev=100000)
    
    with open("output/xstarls.txt",'w') as f:
        f.write("%11.8f\n" % popt[0])
        f.write("%11.8f\n" % popt[1])
        f.write("%11.8f\n" % popt[2])
        f.write("%11.8f" % popt[3])
    
else:
    df_file = "output/original.txt"
    df = pd.read_csv(df_file,header=None, sep=" ")
    popt, pcov = curve_fit(func_original,df[0].values,df[1].values,method='dogbox',maxfev=100000)
