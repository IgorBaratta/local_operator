import pandas as pd
import seaborn
import matplotlib.pyplot as plt
import numpy

dt = pd.read_csv("N1curl.txt")
seaborn.set_context("talk")
rank = 1

unique_flags = dt["flags"].unique()
print(unique_flags)
dt = dt.query(f"flags!='{unique_flags[1]}'").copy()
dt = dt.query(f"flags!='{unique_flags[3]}'").copy()
dt = dt.query(f"flags!='{unique_flags[4]}'").copy()
unique_flags = dt["flags"].unique()
print(unique_flags)

compilers = dt["compiler"].unique()
dt = dt.query(f"rank=={rank}").copy()
dt = dt[dt["compiler"] == compilers[0]]
form_compilers = dt["fcomp"].unique()
dt = dt[dt["fcomp"] == form_compilers[0]]

dt["Time[s]"] = dt["time"]
dt["form compiler"] = dt["fcomp"]

dt["flags"][dt["flags"] == unique_flags[0]] = "AVX512"
dt["flags"][dt["flags"] == unique_flags[1]] = "Scalar"


degrees = dt["degree"].unique()
dt["speed up"] = 0


dt1 = dt.query(f"flags=='AVX512'").copy()
dt2 = dt.query(f"flags=='Scalar'").copy()


avg = numpy.zeros(5, dtype=numpy.float64)
for i in degrees:
    dt2 = dt.query(f"degree=={i}").copy()
    avg[i-1] = dt2["time"].max()

print(avg)

speedup = numpy.zeros(dt1.shape[0])
for i in range(dt1.shape[0]):
    scalar = avg[dt1["degree"].array[i] - 1]
    speedup[i] = scalar/dt1["time"].array[i].min()

dt1["speed up"] = speedup

g = seaborn.catplot(x="degree", y="speed up",
                    kind="bar", data=dt1)
g.set(ylim=(0, 8))
plt.show()
