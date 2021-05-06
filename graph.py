import seaborn
import pandas as pd
import matplotlib.pyplot as plt
import sys

if len(sys.argv) < 2:
    raise RuntimeError("Use with data filename Form.txt")
data_filename = sys.argv[1]

if len(sys.argv) > 2:
    degree = sys.argv[-1]
else:
    degree = 1

dt = pd.read_csv(data_filename)

problem = dt['problem'][0]
assert all(dt['problem'] == problem)
print(problem, degree)

seaborn.set(style="ticks")
seaborn.set_style("darkgrid")

g = seaborn.catplot(x="compiler", y="time", hue="method", col="flags", kind="bar", data=dt)
# g.set(yticks=list(range(5)), ylim=(0, 6))
g.set_titles(problem + " " + str(degree))
plt.show()
