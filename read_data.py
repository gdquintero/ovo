import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

def model1(x,t):
    a = x[0]
    b = x[1]
    c = x[2]
    ebt = np.exp(-1.0 * b * t)

    res = (a / b) * t * ebt + (1.0 / b) * ((a / b) - c) * (ebt - 1.0) - c * t
    return 1.0 - np.exp(res)

def model2(x,t):
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
x[1] = float(xdata[2])
x[2] = float(xdata[1])

tmin = df["age"].values[0]
tmax = df["age"].values[-1]

t = np.linspace(tmin,tmax,1000)

plt.plot(df["age"],df["ratio"],"o")
plt.plot(t,model1(x,t))
plt.show()
plt.close()

fig, ax = plt.subplots()
ax = sns.boxplot(x=df["ratio"],ax=ax)

# plt.show()