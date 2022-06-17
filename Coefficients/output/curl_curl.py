import pandas as pd
import seaborn
import matplotlib.pyplot as plt
import numpy

dt = pd.read_csv("curl_curl.txt")
# seaborn.set(font_scale=1.5)
seaborn.set_context("talk")
seaborn.set(font_scale = 2)
seaborn.set_style("whitegrid")

dt["Time (s)"] =  dt["time"]
dt["Form Compiler"] =  dt["form_compiler"]

g = seaborn.barplot(x="Architecture", y="Time (s)", hue="Form Compiler", data=dt)

plt.grid("both")
plt.show()
