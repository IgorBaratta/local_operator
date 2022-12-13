
from ffcx.codegeneration.backend import FFCXBackend
from ffcx.analysis import analyze_ufl_objects
from ffcx.ir.representation import compute_ir
from ffcx.codegeneration.integrals import IntegralGenerator
from ffcx.element_interface import create_element
from ffcx.codegeneration.C.format_lines import format_indented_lines
from ffcx.options import get_options
import ufl
import typing
from importlib import reload
from compiler import problem


_arguments = """({scalar_type}* restrict A,
                   const {scalar_type}* restrict w,
                   const {scalar_type}* restrict c,
                   const {geom_type}* restrict coordinate_dofs,
                   const int* restrict entity_local_index,
                   const uint8_t* restrict quadrature_permutation)\n"""


_headers = """
#include <cmath>
#include <cstdint>
#include <complex.h>
#define restrict __restrict__

constexpr int dim = {dim};
constexpr int global_size = {global_size};
constexpr int kernel_rank = {rank};
constexpr int num_nodes = {num_nodes};
constexpr int batch_size = {batch_size};
constexpr int num_coefficients = {num_coefficients};


using scalar_type={scalar_type};
using geom_type={geom_type};
using namespace std;
"""


_headers_batched = """
#include <cmath>
#include <cstdint>

#define restrict __restrict__
#define USE_VECTOR_EXTENSIONS


constexpr int dim = {dim};
constexpr int global_size = {global_size};
constexpr int kernel_rank = {rank};
constexpr int num_nodes = {num_nodes};
constexpr int num_coefficients = {num_coefficients};

#if defined(__clang__)
    typedef {scalar_type} {scalar_type}{batch_size} __attribute__((ext_vector_type({batch_size})));
#elif defined(__GNUC__) || defined(__GNUG__)
    typedef {scalar_type} {scalar_type}{batch_size} __attribute__((vector_size({batch_size} * sizeof({scalar_type}))));
#else
    #error "Compiler not supported"
#endif

#if defined(__clang__)
    typedef {geom_type} {geom_type}{batch_size} __attribute__((ext_vector_type({batch_size})));
#elif defined(__GNUC__) || defined(__GNUG__)
    typedef {geom_type} {geom_type}{batch_size} __attribute__((vector_size({batch_size} * sizeof({geom_type}))));
#else
    #error "Compiler not supported"
#endif

using scalar_type = {scalar_type}{batch_size};
using geom_type = {geom_type}{batch_size};

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
        parameters = get_options()

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

    scalar_type = parameters["scalar_type"]
    geom_type = scalar_type.replace(' _Complex', '')
    batch_size = parameters["batch_size"]

    if batch_size and batch_size > 1:
        geom_type += str(batch_size)
        scalar_type += str(batch_size)

    settings = {"scalar_type": scalar_type, "geom_type": geom_type}
    arguments = _arguments.format(**settings)
    signature = "inline void " + name + arguments
    body = compute_integral_body(integral_ir, backend)
    code = signature + " {\n" + body + "\n}\n"

    return code


def generate_code(action, scalar_type, global_size, batch_size):
    reload(problem)

    batch_size = batch_size if batch_size else 1
    parameters = get_options()
    parameters["scalar_type"] = scalar_type
    parameters["batch_size"] = batch_size

    if action:
        code = compile_form(problem.L, "kernel", parameters)
        num_coefficients = analyze_ufl_objects(
            [problem.L], parameters).form_data[0].num_coefficients
        rank = 1
    else:
        code = compile_form(problem.a, "kernel", parameters)
        num_coefficients = analyze_ufl_objects(
            [problem.a], parameters).form_data[0].num_coefficients
        rank = 2

    element = create_element(problem.element)
    num_nodes = element.cell().num_vertices()
    geom_type = scalar_type.replace(' _Complex', '')

    if batch_size > 1:
        headers = _headers_batched.format(dim=element.dim, global_size=global_size,
                                          scalar_type=scalar_type, rank=rank, geom_type=geom_type,
                                          batch_size=batch_size, num_nodes=num_nodes,
                                          num_coefficients=num_coefficients)
    else:
        headers = _headers.format(dim=element.dim, global_size=global_size,
                                  scalar_type=scalar_type, rank=rank, geom_type=geom_type,
                                  batch_size=batch_size, num_nodes=num_nodes, num_coefficients=num_coefficients)

    with open("compiler/problem.hpp", "w") as file:
        file.write(headers)
        file.write(code)