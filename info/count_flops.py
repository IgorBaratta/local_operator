from typing import Optional

import ffcx.options
import ufl
from ffcx.analysis import analyze_ufl_objects
from ffcx.codegeneration.backend import FFCXBackend
from ffcx.codegeneration.integrals import IntegralGenerator
from ffcx.ir.representation import compute_ir
import numpy


def count_flops(form: ufl.Form, options: Optional[dict] = {}):
    """Return a list with the number of flops for each kernel in the Form."""
    options = ffcx.options.get_options(options)
    assert isinstance(form, ufl.Form)
    analysis = analyze_ufl_objects([form], options)
    ir = compute_ir(analysis, {}, "flops", options, False)

    flops = []
    _bytes = []

    for integral_ir in ir.integrals:
        # Create FFCx C backend
        backend = FFCXBackend(integral_ir, options)
        # Configure kernel generator
        ig = IntegralGenerator(integral_ir, backend)
        # Generate code ast for the tabulate_tensor body
        ast = ig.generate()
        _sum = 0

        for statement in ast.statements:
            if isinstance(statement, ffcx.codegeneration.C.cnodes.ArrayDecl):
                _sum += numpy.prod(statement.sizes)
            if isinstance(statement, ffcx.codegeneration.C.cnodes.Scope):
                for sub_statements in statement.body.statements:
                    if isinstance(sub_statements, ffcx.codegeneration.C.cnodes.ArrayDecl):
                        _sum += numpy.prod(sub_statements.sizes)
        flops.append(ast.flops())
        _bytes.append(_sum)
        
    return flops[0], _bytes[0]
