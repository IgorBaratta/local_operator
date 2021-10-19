from ufl import *

degree = 4

element = FiniteElement("Lagrange", tetrahedron, degree)
mesh = Mesh(VectorElement("Lagrange", tetrahedron, 1))

V = FunctionSpace(mesh, element)
u = TrialFunction(V)
v = TestFunction(V)

W = FunctionSpace(mesh, FiniteElement("Lagrange", tetrahedron, 1))
k = Coefficient(W)

a = k*inner(grad(u), grad(v))*dx