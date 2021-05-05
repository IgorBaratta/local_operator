import seaborn
import pandas as pd
import matplotlib.pyplot as plt
import numpy
import sys

if len(sys.argv) > 1:
    degree = sys.argv[1]
else:
    degree = 0

dt = pd.read_csv("Lagrange" + degree)

seaborn.set(style="ticks")
seaborn.set_style("darkgrid")

g = seaborn.catplot(x="compiler", y="time", hue="method", col="flags", kind="bar", data=dt)
# g.set(yticks=list(range(5)), ylim=(0, 6))
g.set_titles("Lagrange" + degree)
plt.show()