import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

def model(x,t):
    a = x[0]
    b = x[1]
    c = x[2]
    ebt = np.exp(-1.0 * b * t)

    res = (a / b) * t * ebt + (1.0 / b) * ((a / b) - c) * (ebt - 1.0) - c * t
    return 1.0 - np.exp(res)

def least_squares(t):
    a = 0.051
    b = 0.33
    c = 0.003
    ebt = np.exp(-1.0 * b * t)

    res = (a / b) * t * ebt + (1.0 / b) * ((a / b) - c) * (ebt - 1.0) - c * t
    return 1.0 - np.exp(res)

def lamb(x,t):
    a = x[0]
    b = x[1]
    c = x[2]
    ebt = np.exp(-1.0 * b * t)

    return (a * t - c) * ebt + c

df = pd.read_excel("zika.xlsx")

with open("output/zika.txt","w") as f:
    for i in range(len(df["age"].values)):
        f.write("%f %f\n" % (df["age"].values[i],df["ratio"].values[i]))

with open("output/xstarovo.txt") as f:
    lines = f.readlines()
    xdata = [line.split()[0] for line in lines]

x = np.empty(3)


x[0] = float(xdata[0])
x[1] = float(xdata[1])
x[2] = float(xdata[2])

tmin = df["age"].values[0]
tmax = df["age"].values[-1]

t = np.linspace(tmin,tmax,1000)

fig, axs = plt.subplots(2,1)

axs[0].plot(df["age"],df["ratio"],"ko")
axs[0].plot(t,least_squares(t),label="LS",color="red")
axs[0].plot(t,model(x,t),label="OVO",color="blue")
axs[0].legend()

axs[1].plot(t,lamb([0.051,0.33,0.003],t),label="LS",color="red")
axs[1].plot(t,lamb(x,t),label="OVO",color="blue")
axs[1].legend()

plt.show()
plt.close()

fig, ax = plt.subplots()
ax = sns.boxplot(data=df,x="ratio",whis=1.0,ax=ax)
# plt.show()