import problem
import tsfc
import loopy
from ffcx.element_interface import create_element


def generate_code(matrix_free: bool):

    if matrix_free:
        form, = tsfc.compile_form(problem.L, coffee=False)
        rank = 1
    else:
        form, = tsfc.compile_form(problem.a, coffee=False)
        rank = 2

    ffcx_element = create_element(problem.element)

    headers = "#include <cmath>"
    headers += "\n#include <cstdint>"
    headers += "\n\n#define restrict __restrict__"
    headers += f"\n\n#define MF {rank}"
    headers += "\n\nusing namespace std;"
    headers += f"\n\nconstexpr int dim = {ffcx_element.dim};"
    headers += f"\n\nconstexpr int kernel_rank = {rank}; \n\n"

    kernel = loopy.generate_code_v2(form.ast).device_code()
    kernel = kernel.replace("void", "static inline void")

    with open("tsfc/problem.hpp", "w") as file:
        file.write(headers)
        file.write(kernel)
