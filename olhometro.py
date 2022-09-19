import numpy as np
import pandas as pd
from sklearn.metrics import mean_squared_error

def func(t,a,b,c):
    ebt = np.exp(-1.0 * b * t)
    res = (a / b) * t * ebt + (1.0 / b) * ((a / b) - c) * (ebt - 1.0) - c * t
    res = 1.0 - np.exp(res)

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

error_ovo = mean_squared_error(df["ratio"].values,func(df["age"].values,*x_ovo))
error_ls = mean_squared_error(df["ratio"].values,func(df["age"].values,*x_ls))

print(error_ovo,error_ls)

y_true = np.empty(len(df["age"]) - 4)
y_pred_ovo = np.empty(len(df["age"]) - 4)
y_pred_ls = np.empty(len(df["age"]) - 4)

y_true[:15] = df["ratio"].values[:15]
y_true[15:] = df["ratio"].values[19:]

y_pred_ovo[:15] = func(df["age"].values[:15],*x_ovo)
y_pred_ovo[15:] = func(df["age"].values[19:],*x_ovo)

y_pred_ls[:15] = func(df["age"].values[:15],*x_ls)
y_pred_ls[15:] = func(df["age"].values[19:],*x_ls)

error_ovo_outliers = mean_squared_error(y_true,y_pred_ovo)
error_ls_outliers = mean_squared_error(y_true,y_pred_ls)

print(error_ovo_outliers,error_ls_outliers)