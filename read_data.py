import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

df = pd.read_excel("zika.xlsx")

with open("output/zika.txt","w") as f:
    for i in range(len(df["age"].values)):
        f.write("%f %f\n" % (df["age"].values[i],df["ratio"].values[i]))


plt.plot(df["age"],df["ratio"],"o")
plt.show()
plt.close()

fig, ax = plt.subplots()
ax = sns.boxplot(x=df["ratio"],ax=ax)

plt.show()