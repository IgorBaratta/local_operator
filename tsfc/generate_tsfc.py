import problem
import tsfc
import loopy
from ffcx.element_interface import create_element
from ffcx.naming import cdtype_to_numpy
from tsfc.parameters import PARAMETERS as tsfc_default_parameters
from tsfc.kernel_interface.firedrake import Kernel
from importlib import reload
import numpy


_headers = """
#include <cmath>
#include <cstdint>
#define restrict __restrict__

constexpr int dim = {dim};
constexpr int global_size = {global_size};
constexpr int kernel_rank = {rank};

using scalar_type={scalar_type};
using namespace std;
"""


def generate_code(action: bool, scalar_type: str, global_size: int):
    reload(problem)

    tsfc_default_parameters["scalar_type"] = cdtype_to_numpy(scalar_type)
    tsfc_default_parameters["scalar_type_c"] = scalar_type

    if action:
        form, = tsfc.compile_form(
            problem.L, coffee=False, parameters=tsfc_default_parameters)
        rank = 1
    else:
        form, = tsfc.compile_form(problem.a, coffee=False)
        rank = 2

    ffcx_element = create_element(problem.element)

    headers = _headers.format(dim=ffcx_element.dim, global_size=global_size,
                              scalar_type=scalar_type, rank=rank)


    if isinstance(form, Kernel):
        kernel = form.ast.gencode()
    else:
        kernel = loopy.generate_code_v2(form.ast).device_code()
        kernel = kernel.replace("void", "static inline void")

    with open("tsfc/problem.hpp", "w") as file:
        file.write(headers)
        file.write(kernel)
