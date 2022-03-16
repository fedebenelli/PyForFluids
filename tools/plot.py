import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("dat2", sep=" ")
df_gpec = pd.read_csv("gpeciso", sep=" ")
df_gpec322 = pd.read_csv("gpeciso322", sep=" ")

axes = df.groupby("x").plot.scatter(x="it", y="p")
axes = df.groupby("x").plot.scatter(x="it", y="ki")
for x in axes.index:
    axes[x].set_title(x)
plt.show()


gerg_bub = df.groupby("x").tail(1)

plt.scatter(gerg_bub["x"], gerg_bub["p"], alpha=0.5)
plt.scatter(gerg_bub["y"], gerg_bub["p"], alpha=0.5)

plt.plot(df_gpec["x"], df_gpec["p"]*1e5)
plt.plot(df_gpec["y"], df_gpec["p"]*1e5)

plt.plot(df_gpec322["x"], df_gpec322["p"]*1e5)
plt.plot(df_gpec322["y"], df_gpec322["p"]*1e5)

plt.ylim(0, 1e7)

plt.show()
