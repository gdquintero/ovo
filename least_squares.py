import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

def func(t,a,b,c):
    ebt = np.exp(-1.0 * b * t)
    res = (a / b) * t * ebt + (1.0 / b) * ((a / b) - c) * (ebt - 1.0) - c * t
    res = 1.0 - np.exp(res)

    return res

outlier = True

if outlier:
    df_file = "zika_outliers.xlsx"
else:
    df_file = "zika.xlsx"

df = pd.read_excel(df_file)

popt, pcov = curve_fit(func,df["age"].values,df["ratio"].values,bounds=(0,np.inf * np.ones(3)))

with open("output/xstarls.txt",'w') as f:
    f.write("%11.8f\n" % popt[0])
    f.write("%11.8f\n" % popt[1])
    f.write("%11.8f" % popt[2])

print(popt)