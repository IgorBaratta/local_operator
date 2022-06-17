import pandas as pd
import seaborn
import matplotlib.pyplot as plt
import numpy

dt = pd.read_csv("N1curl.txt")

degree = 4
rank = 2

compilers = dt["compiler"].unique()
unique_flags = dt["flags"].unique()
dt = dt.query(f"degree=={degree}").copy()
dt = dt.query(f"rank=={rank}").copy()
dt = dt[dt["flags"] == unique_flags[0]]
dt = dt[dt["compiler"] == compilers[0]]


dt["Time[s]"] = dt["time"]
g = seaborn.catplot(x="compiler", y="Time[s]", hue="fcomp",
                    kind="bar", data=dt)
g.fig.subplots_adjust(top=.85)
g.fig.suptitle(
    f"k*inner(curl(u), curl(v))*dx \n degree = {degree} $\,$ rank = {rank}")
g._legend.remove()
plt.tight_layout()
plt.legend(loc="best")
plt.show()


# degrees = dt["degree"].unique()
# degree = numpy.arange(1, 6, dtype=int)
# form_compilers = ["' ffcx'", "' tsfc'"]#, "' ffc'"]

# fig, ax = plt.subplots()

# for form_compiler in form_compilers:
#     aux = dt.groupby(['fcomp', 'degree']).mean()
#     time = aux.query(f"fcomp == {form_compiler}").time.to_numpy()
#     dofs = (degree + 1)*(degree + 2) * (degree + 3)/6 * 500000
#     metric = dofs/time
#     ax.plot(degree, metric, label=form_compiler)


# ax.legend()
# ax.set_yscale('log')
# ax.set(xlabel='Polynomial degree', ylabel='Dofs/s')
# ax.set(title="Elasticity - Rank 1")
# ax.grid()
# plt.show()
