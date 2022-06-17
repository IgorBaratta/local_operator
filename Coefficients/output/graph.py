import pandas as pd
import seaborn
import matplotlib.pyplot as plt
import numpy

dt = pd.read_csv("blas.txt")
seaborn.set_context("talk")

cells = 10000
dt["Dofs/s"] = dt["dofs"]/dt["time"]
dt["Cells/s"] = 10000/dt["time"]
dt["Polynomial Degree"] = dt["degree"]
dt["Method"] = dt["method"]
dt["Cell Type"] = dt["celltype"]


dt = dt.query("celltype == ' hexahedron'")
dt = dt.query("Method != ' Blis'").copy()

g = seaborn.lineplot(x="Polynomial Degree", y="Cells/s", hue="Method",
                     data=dt, style="Method", markers=True, linewidth=4)
g.set_yscale("log")
g.set_xscale("log")
g.legend(bbox_to_anchor=(0.7, 1.05))


C = 2e-7
degrees = numpy.arange(3, 11)
cost = numpy.zeros(degrees.size)
for i in range(degrees.size):
    cost[i] = cells/(C*(degrees[i]+1)**(4))

plt.plot(degrees, cost, "k--", linewidth=4, label=r"$p^{-4}$")


C = 1e-6
for i in range(degrees.size):
    cost[i] = cells/(C*(degrees[i]+1)**(6))
plt.plot(degrees, cost, "k-.", label=r"$p^{-6}$", linewidth=4)

plt.legend()
plt.grid(True, which="both", ls=":")
plt.show()
