from ufl import *

element = VectorElement("Lagrange", tetrahedron, 3)
mesh = Mesh(VectorElement("Lagrange", tetrahedron, 1))

V = FunctionSpace(mesh, element)
u = TrialFunction(V)
v = TestFunction(V)

a = inner(u, v)*dx

un = Coefficient(V)
L = action(a, un)