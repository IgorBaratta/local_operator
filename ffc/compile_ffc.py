import problem
import os
import ffc
from ffc.fiatinterface import create_element
from importlib import reload


def generate_code(matrix_free: bool):
    reload(problem)
    if matrix_free:
        rank = 1
    else:
        rank = 2

    os.environ["UFC_INCLUDE_DIR"] = ffc.get_include_path()

    if os.system(f"cd ffc && ffc problem.ufl") != 0:
        raise RuntimeError("ffc failed")

    ffc_element = create_element(problem.element)
    ffc_element.dim = ffc_element.space_dimension()

    headers = "#include <cmath>"
    headers += "\n#include <cstdint>"
    headers += "\n\n#define restrict __restrict__"
    headers += f"\n\n#define MF {rank}"
    headers += "\n\nusing namespace std;"
    headers += f"\n\nconstexpr int dim = {ffc_element.dim};"
    headers += f"\n\nconstexpr int kernel_rank = {rank}; \n\n"

    with open("ffc/data.hpp", "w") as file:
        file.write(headers)
