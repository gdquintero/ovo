import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt

def detect_outlier(data):
    
    threshold = 2
    mean = np.mean(data)
    std = np.std(data)
    
    for y in data:
        z_score= (y - mean)/std 
        if np.abs(z_score) > threshold:
            outliers.append(y)
    return outliers

outlier = True

if outlier:
    df_file = "zika_outliers.xlsx"
    txt_file = "output/zika_outliers.txt"
else:
    df_file = "zika.xlsx"
    txt_file = "output/zika.txt" 

df = pd.read_excel(df_file)
outliers = []

plt.plot(df["age"],df["ratio"],"ko")
plt.show()
plt.close()

with open(txt_file,"w") as f:
    for i in range(len(df["age"].values)):
        f.write("%f %f\n" % (df["age"].values[i],df["ratio"].values[i]))

fig, ax = plt.subplots()
ax = sns.boxplot(data=df,x="ratio",whis=1.0,ax=ax)
# plt.show()

outliers = detect_outlier(df["ratio"])
