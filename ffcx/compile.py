
from ffcx.codegeneration.backend import FFCXBackend
from ffcx.analysis import analyze_ufl_objects
from ffcx.ir.representation import compute_ir
from ffcx.codegeneration.integrals import IntegralGenerator
from ffcx.element_interface import create_element
from ffcx.codegeneration.C.format_lines import format_indented_lines
from ffcx import get_parameters
import basix
import ufl
import typing
import problem


def compute_integral_body(ir, backend):

    # Configure kernel generator
    ig = IntegralGenerator(ir, backend)

    # Generate code ast for the tabulate_tensor body
    parts = ig.generate()

    # Format code as string
    body = format_indented_lines(parts.cs_format(ir.precision), 1)

    return body


def compile_form(form: ufl.Form, name: str,
                 parameters: typing.Dict = None,
                 visualise: bool = False):

    if parameters is None:
        parameters = get_parameters()

    # Stage 1: analysis
    analysis = analyze_ufl_objects([form], parameters)

    # Stage 2: intermediate representation
    ir = compute_ir(analysis, {}, " ", parameters, visualise)

    if len(ir.integrals) > 1:
        raise RuntimeError(
            "This function is meant to compile one integral type a time.")

    # Stage 3: code generation
    integral_ir = ir.integrals[0]
    backend = FFCXBackend(integral_ir, parameters)

    # FIXME: Get signature from backend (c++/template, c, cuda, etc)
    arguments = """(double* restrict A,
                    const double* restrict w,
                    const double* restrict c,
                    const double* restrict coordinate_dofs,
                    const int* restrict entity_local_index,
                    const uint8_t* restrict quadrature_permutation) {\n"""

    signature = "inline void " + name + arguments
    body = compute_integral_body(integral_ir, backend)
    code = signature + body + "\n}\n"

    return code


def generate_code(matrix_free):
    if matrix_free:
        code = compile_form(problem.L, "kernel")
        rank = 1
    else:
        code = compile_form(problem.a, "kernel")
        rank = 2

    ffcx_element = create_element(problem.element)

    headers = "#include <cmath>"
    headers += "\n#include <cstdint>"
    headers += "\n\n#define restrict __restrict__"
    headers += "\n\nusing namespace std;"
    headers += f"\n\nconstexpr int dim = {ffcx_element.dim};"
    headers += f"\n\nconstexpr int kernel_rank = {rank}; \n\n"

    with open("ffcx/problem.hpp", "w") as file:
        file.write(headers)
        file.write(code)
