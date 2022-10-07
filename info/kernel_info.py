import os
import platform
import yaml
from string import Template
from subprocess import Popen, PIPE
from importlib import reload, import_module
from count_flops import count_flops
import numpy
from matplotlib import pyplot as plt
from ffcx.element_interface import create_element


problem_ = "Laplacian"
degrees = numpy.arange(1, 15)
cell_type = "hexahedron"

for degree in degrees:
    with open("../forms/" + problem_ + ".ufl", 'r') as f:
        src = Template(f.read())
        d = {'degree': str(degree), 'vdegree': str(
            degree + 1), "cell": cell_type}
        result = src.substitute(d)

    with open(f"problem{degree}.py", "w") as f2:
        f2.writelines(result)


flops_per_dof = numpy.zeros_like(degrees)
bytes_per_dof = numpy.zeros_like(degrees)
cache_per_dof = numpy.zeros_like(degrees)
dofs_list = numpy.zeros_like(degrees)
for i, degree in enumerate(degrees):
    module = import_module(f"problem{degree}")
    form = module.L
    element = create_element(module.element)
    dofs = element.dim
    dofs_list[i] = dofs

    cell = element.cell()
    geometry_size = cell.num_vertices() * 3

    coeff_size = 0
    coefficients = form.coefficients()
    for coeff in coefficients:
        el = create_element(coeff.ufl_element())
        coeff_size += el.dim

    flops_per_dof[i], cache_per_dof[i] = count_flops(form)
    bytes_per_dof[i] = ((coeff_size + 2*dofs + geometry_size) * 8)/dofs
    cache_per_dof[i] = (cache_per_dof[i] * 8)/dofs
    flops_per_dof[i] = flops_per_dof[i]/dofs


FLOPS = numpy.full_like(degrees, 5e12, dtype=float)
BW = numpy.full_like(degrees, 300e9, dtype=float)
BW1 = numpy.full_like(degrees, 2.0e4*1e9, dtype=float)
BW2 = numpy.full_like(degrees, 0.5e4*1e9, dtype=float)
BW3 = numpy.full_like(degrees, 2.3e3*1e9, dtype=float)


cache_rate = numpy.zeros_like(BW)
level = numpy.zeros_like(BW)

for i, num_dofs in enumerate(dofs_list):
    size = num_dofs * (cache_per_dof[i] + bytes_per_dof[i])
    ai_cache = numpy.divide(flops_per_dof[i], size/num_dofs)
    if size*1.5 < 32e3:
        cache_rate[i] = ai_cache * BW1[i]
        level[i] = 1
    elif size*1.5 < 1280e3:
        cache_rate[i] = ai_cache * BW2[i]
        level[i] = 2
    else:
        cache_rate[i] = ai_cache * BW3[i]
        level[i] = 3


max_throughput_cache = cache_rate/flops_per_dof
ai_std = flops_per_dof/bytes_per_dof
max_flops = numpy.minimum(FLOPS, ai_std*BW)
max_throughput = max_flops/flops_per_dof

markers = ["o", "v", "s", "X"]
plt.plot(degrees, max_throughput, label="standard model",
         marker=markers[1], linestyle='dashed')

max_throughput_cache = numpy.minimum(max_throughput, max_throughput_cache)

plt.plot(degrees, max_throughput_cache, label="cache-aware model",
         marker=markers[1], linestyle='dashed')


plt.ylabel("max throughput(dofs/s)")
plt.legend()
plt.xlabel(r"polynomial degree $P$")
plt.grid(True, which="both")
plt.yscale("log")
plt.show()


print(level)
