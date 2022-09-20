import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import mean_squared_error

def func(t,a,b,c):
    ebt = np.exp(-1.0 * b * t)
    res = (a / b) * t * ebt + (1.0 / b) * ((a / b) - c) * (ebt - 1.0) - c * t
    res = 1.0 - np.exp(res)

    return res

def model(t,a,b,c):
    ebt = np.exp(-1.0 * b * t)
    res = (a * t - c) * ebt + c

    return res

with open("output/xstarovo.txt") as f:
    lines = f.readlines()
    xdata = [line.split()[0] for line in lines]

x_ovo = np.empty(3)
x_ls = np.empty(3)

for i in range(3):
    x_ovo[i] = float(xdata[i])

with open("output/xstarls.txt") as f:
    lines = f.readlines()
    xdata = [line.split()[0] for line in lines]

for i in range(3):
    x_ls[i] = float(xdata[i])

outlier = True

if outlier:
    df_file = "zika_outliers.xlsx"
else:
    df_file = "zika.xlsx"

df = pd.read_excel(df_file)

y_true = np.empty(len(df["age"]) - 4)
y_pred_ovo = np.empty(len(df["age"]) - 4)
y_pred_ls = np.empty(len(df["age"]) - 4)

y_true[:15] = df["ratio"].values[:15]
y_true[15:] = df["ratio"].values[19:]

y_pred_ovo[:15] = func(df["age"].values[:15],*x_ovo)
y_pred_ovo[15:] = func(df["age"].values[19:],*x_ovo)

y_pred_ls[:15] = func(df["age"].values[:15],*x_ls)
y_pred_ls[15:] = func(df["age"].values[19:],*x_ls)

error_ovo = mean_squared_error(y_true,y_pred_ovo)
error_ls = mean_squared_error(y_true,y_pred_ls)

fig, ax = plt.subplots()

ax.plot(df["age"].values,df["ratio"],"ko")
lines = []
lines = ax.plot(df["age"].values,func(df["age"].values,*x_ovo),"b")
lines += ax.plot(df["age"].values,func(df["age"].values,*x_ls),"r")
ax.legend(lines[:],['OVO', 'Least Squares'],loc='upper right', frameon=False)

textstr = '\n'.join((
    "Error OVO: %.6f" % error_ovo,
    "Error LS: %.6f" % error_ls)
)

plt.text(0.5, 0.96, textstr,
         ha="left", va="center",
         bbox=dict(boxstyle="round",
                   ec="#D5D2D2",
                   fc="#EEEEEE",
                   )
         )
         
plt.show()
plt.close()

t = np.linspace(0,35,1000)

plt.plot(t,model(t,*x_ovo),"b",label="OVO")
plt.plot(t,model(t,*x_ls),"r",label="LS")
plt.legend()
plt.show()