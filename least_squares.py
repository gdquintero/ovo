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

outlier = True
model = 1

if model == 0:
    if outlier:
        df_file = "zika_outliers.xlsx"
    else:
        df_file = "zika.xlsx"

    df = pd.read_excel(df_file)
    popt, pcov = curve_fit(func_zika,df["age"].values,df["ratio"].values,bounds=(0,np.inf * np.ones(3)))

    with open("output/xstarls.txt",'w') as f:
        f.write("%11.8f\n" % popt[0])
        f.write("%11.8f\n" % popt[1])
        f.write("%11.8f" % popt[2])
else:
    df_file = "output/data_wave.txt"
    df = pd.read_csv(df_file,header=None, sep=" ")
    popt, pcov = curve_fit(func_wave,df[0].values,df[1].values,method='dogbox',maxfev=100000)
    
    with open("output/xstarls.txt",'w') as f:
        f.write("%11.8f\n" % popt[0])
        f.write("%11.8f\n" % popt[1])
        f.write("%11.8f\n" % popt[2])
        f.write("%11.8f" % popt[3])


t = np.linspace(0,1.5 * np.pi,1000)

plt.plot(df[0].values,df[1].values,"o")
plt.plot(t,func_wave(t,*popt))
plt.show()
