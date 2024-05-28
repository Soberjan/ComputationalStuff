from dolfin import *
import mpi4py.MPI

comm = mpi4py.MPI.COMM_WORLD
mesh = UnitDiscMesh.create(comm, 10, 2, 2)

# Build function space with Lagrange multiplier
P1 = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
R = FiniteElement("Real", mesh.ufl_cell(), 0)
W = FunctionSpace(mesh, P1 * R)

# Define variational problem
(u, c) = TrialFunction(W)
(v, d) = TestFunctions(W)
eps = 0.1
t = 0
f1 = Expression('x[0] < -0.5 + eps && x[0] > -0.5 - eps && x[1] < -0.5 + eps && x[1] > -0.5 - eps ? sin(t) : 0', eps=eps, t=t, degree=1)
f2 = Expression('x[0] < 0.5 + eps && x[0] > 0.5 - eps && x[1] < 0.5 + eps && x[1] > 0.5 - eps ? cos(t) : 0', eps=eps, t=t, degree=1)
g = Expression("0", degree=2)
a = (inner(grad(u), grad(v)) + c * v + u * d) * dx
L = (f1 + f2) * v * dx

# Compute solution
w = Function(W)
solve(a == L, w)
(u, c) = w.split()

# Save solution in VTK format
file = File("neumann_poisson.pvd")
file << u
