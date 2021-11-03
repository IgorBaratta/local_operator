from ufl import *

P2 = VectorElement("Lagrange", tetrahedron, 4)
P1 = FiniteElement("Lagrange", tetrahedron, 3)
TH = P2 * P1


(u, p) = TrialFunctions(TH)
(v, q) = TestFunctions(TH)

a = (inner(grad(u), grad(v)) - div(v)*p + div(u)*q)*dx


(un, pn) = Coefficients(TH)
L = (inner(grad(un), grad(v)) - div(v)*pn + div(un)*q)*dx