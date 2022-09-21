import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import mean_squared_error

def func(t,a,b,c):
    ebt = np.exp(-1.0 * b * t)
    res = (a / b) * t * ebt + (1.0 / b) * ((a / b) - c) * (ebt - 1.0) - c * t
    res = 1.0 - np.exp(res)

    return res

def func2(t,a,b,c):
    ebt = np.exp(-1.0 * b * t)
    res = (a * t - c) * ebt + c

    return res

def wave(t,a,b,c,d):
    return a * np.exp(b * np.sin(c * t)) + d

with open("output/xstarovo.txt") as f:
    lines = f.readlines()
    xdata = [line.split()[0] for line in lines]

model = 1
outlier = True

if model == 0:
    n = 3
    outliers = 4
    t = np.linspace(0,35,1000)
else:
    n = 4
    outliers = 3
    t = np.linspace(0,1.5 * np.pi,1000)

x_ovo = np.empty(n)
x_ls = np.empty(n)

for i in range(n):
    x_ovo[i] = float(xdata[i])

with open("output/xstarls.txt") as f:
    lines = f.readlines()
    xdata = [line.split()[0] for line in lines]

for i in range(n):
    x_ls[i] = float(xdata[i])

if model == 0:
    if outlier:
        df_file = "output/zika_outliers.txt"
    else:
        df_file = "output/zika.txt"

    df = pd.read_csv(df_file,header=None, sep=" ")

    y_true = np.empty(len(df) - outliers)
    y_pred_ovo = np.empty(len(df) - outliers)
    y_pred_ls = np.empty(len(df) - outliers)

    y_true[:15] = df[1].values[:15]
    y_true[15:] = df[1].values[19:]

    y_pred_ovo[:15] = func(df[0].values[:15],*x_ovo)
    y_pred_ovo[15:] = func(df[0].values[19:],*x_ovo)

    y_pred_ls[:15] = func(df[0].values[:15],*x_ls)
    y_pred_ls[15:] = func(df[0].values[19:],*x_ls)
else:
    df_file = "output/wave.txt"
    df = pd.read_csv(df_file,header=None, sep=" ")
    
    y_true = np.empty(len(df) - outliers)
    y_pred_ovo = np.empty(len(df) - outliers)
    y_pred_ls = np.empty(len(df) - outliers)

    with open("output/outlier_indices_wave.txt") as f:
        lines = f.readlines()
        outliers = [line.split()[0] for line in lines]
        
    for i in range(3):
        outliers[i] = int(outliers[i])

    j = 0

    for i in range(df[0].size):
        if not i in outliers:
            y_true[j] = df[1].values[i]
            y_pred_ovo[j] = wave(df[0].values[i],*x_ovo)
            y_pred_ls[j] = wave(df[0].values[i],*x_ls)
            j += 1

error_ovo = mean_squared_error(y_true,y_pred_ovo)
error_ls = mean_squared_error(y_true,y_pred_ls)

fig, ax = plt.subplots()

ax.plot(df[0].values,df[1].values,"ko")
lines = []
if model == 0:
    lines = ax.plot(t,func(t,*x_ovo),"b")
    lines += ax.plot(t,func(t,*x_ls),"r")
else:
    lines = ax.plot(t,wave(t,*x_ovo),"b")
    lines += ax.plot(t,wave(t,*x_ls),"r")

ax.legend(lines[:],['OVO', 'Least Squares'],loc='lower right', frameon=False)

# textstr = '\n'.join((
#     "Error OVO: %.6f" % error_ovo,
#     "Error LS: %.6f" % error_ls)
# )

# plt.text(0.5, 0.96, textstr,
#          ha="left", va="center",
#          bbox=dict(boxstyle="round",
#                    ec="#D5D2D2",
#                    fc="#EEEEEE",
#                    )
#          )
         
plt.show()
# plt.close()

# plt.plot(t,func2(t,*x_ovo),"b",label="OVO")
# plt.plot(t,func2(t,*x_ls),"r",label="LS")
# plt.legend()
# plt.show()