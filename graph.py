import seaborn
import pandas as pd
import matplotlib.pyplot as plt
import numpy


dt = pd.read_csv("results_archer")
matrix_name = "Mass"
# seaborn.set_theme(style="darkgrid")
# seaborn.set_theme(style="ticks")

seaborn.set(style="ticks")
seaborn.set_style("darkgrid")

dt1 = dt[dt["method"] == "ffcx"].copy()
dt2 = dt[dt["method"] == "fused"].copy()

dt1["speedup"] = dt1.time/dt2.time.array


g = seaborn.catplot(x="degree", y="speedup", hue="compiler", col="flags", kind="bar", data=dt1)
g.set(yticks=list(range(5)), ylim=(0, 5))
plt.show()
