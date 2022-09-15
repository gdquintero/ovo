import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

df = pd.read_excel("zika.xlsx")

with open("output/zika.txt","w") as f:
    for i in range(len(df["age"].values)):
        f.write("%f %f\n" % (df["age"].values[i],df["ratio"].values[i]))

fig, ax = plt.subplots()
ax = sns.boxplot(data=df,x="ratio",whis=1.0,ax=ax)
plt.show()