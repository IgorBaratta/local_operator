
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
from importlib import reload


_arguments = """({scalar_type}{batch_size}* restrict A,
                const {scalar_type}{batch_size}* restrict w,
                const {scalar_type}{batch_size}* restrict c,
                const {scalar_type}{batch_size}* restrict coordinate_dofs,
                const int* restrict entity_local_index,
                const uint8_t* restrict quadrature_permutation)"""

_headers = """
#include <cmath>
#include <cstdint>
#define restrict __restrict__

constexpr int dim = {dim};
constexpr int global_size = {global_size};
constexpr int kernel_rank = {rank};

#if defined(__clang__)
    typedef {scalar_type} {scalar_type}{batch_size} __attribute__((ext_vector_type({batch_size})));
#elif defined(__GNUC__) || defined(__GNUG__)
    typedef {scalar_type} {scalar_type}{batch_size} __attribute__((vector_size({batch_size} * sizeof({scalar_type}))));
#else
    #error "Compiler not supported"
#endif

using scalar_type = {scalar_type}{batch_size};
constexpr int batch_size = {batch_size};
using namespace std;
"""


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

    batch_size = parameters["batch_size"]
    scalar_type = parameters["scalar_type"]

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

    arguments = _arguments.format(
        scalar_type=scalar_type, batch_size=batch_size)

    signature = "inline void " + name + arguments + "{"
    body = compute_integral_body(integral_ir, backend)
    code = signature + body + "\n}\n"

    return code


def generate_code(action, scalar_type, global_size, batch_size):
    reload(problem)

    parameters = get_parameters()
    parameters["scalar_type"] = scalar_type
    parameters["batch_size"] = batch_size


    if action:
        code = compile_form(problem.L, "kernel", parameters)
        rank = 1
    else:
        code = compile_form(problem.a, "kernel", parameters)
        rank = 2

    ffcx_element = create_element(problem.element)

    headers = _headers.format(dim=ffcx_element.dim, global_size=global_size,
                              scalar_type=scalar_type, rank=rank, batch_size=batch_size)

    with open("cross/problem.hpp", "w") as file:
        file.write(headers)
        file.write(code)
