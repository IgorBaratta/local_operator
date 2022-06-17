# import FIAT, finat
# from ufl import *

# def gauss_lobatto_legendre_line_rule(degree):
#     fiat_make_rule = FIAT.quadrature.GaussLobattoLegendreQuadratureLineRule
#     fiat_rule = fiat_make_rule(FIAT.ufc_simplex(1), degree + 1)
#     finat_ps = finat.point_set.GaussLobattoLegendrePointSet
#     points = finat_ps(fiat_rule.get_points())
#     weights = fiat_rule.get_weights()
#     return finat.quadrature.QuadratureRule(points, weights)

# def gauss_lobatto_legendre_cube_rule(dimension, degree):
#     make_tensor_rule = finat.quadrature.TensorProductQuadratureRule
#     result = gauss_lobatto_legendre_line_rule(degree)
#     for _ in range(1, dimension):
#         line_rule = gauss_lobatto_legendre_line_rule(degree)
#         result = make_tensor_rule([result, line_rule])
#     return result



# p = 5
# element = FiniteElement("Lagrange", hexahedron, 5)
# mesh = Mesh(VectorElement("Lagrange", hexahedron, 1))
# V = FunctionSpace(mesh, element)
# u = Coefficient(V)
# v = TestFunction(V)
# a = dot(u, v) *dx  # Laplace operator


# from tsfc import compile_form
# import loopy

# kernel_spectral, = compile_form(a)


# gll_quadrature_rule = gauss_lobatto_legendre_cube_rule(dimension=3, degree=p)
# a_gll = dot(grad(u), grad(v)) *dx(rule=gll_quadrature_rule)
# kernel_gll, = compile_form(a_gll, coffee=False)

# kernel = loopy.generate_code_v2(kernel_gll.ast).device_code()
# kernel = kernel.replace("void", "static inline void")
# print(kernel)


from mpi4py import MPI
import dolfinx.generation
from dolfinx.cpp.mesh import CellType

mesh = dolfinx.generation.BoxMesh(MPI.COMM_WORLD,[[0.0,0.0,0.0], [200, 200, 200]], [96, 96, 96], CellType.hexahedron)