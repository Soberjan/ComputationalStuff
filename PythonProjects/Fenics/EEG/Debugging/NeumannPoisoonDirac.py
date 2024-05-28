from dolfin import *
import mpi4py.MPI

comm = mpi4py.MPI.COMM_WORLD
mesh = UnitDiscMesh.create(comm, 10, 2, 2)

# Build function space with Lagrange multiplier
P1 = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
R = FiniteElement("Real", mesh.ufl_cell(), 0)
W = FunctionSpace(mesh, P1 * R)

# Define variational problem
(u, c) = TrialFunctions(W)
(v, d) = TestFunctions(W)

a = (inner(grad(u), grad(v)) + c * v + u * d) * dx
L = Constant(0) * v * dx
A, b = assemble_system(a, L)

delta = PointSource(W, Point(0, 0), 1)
delta.apply(b)
delta = PointSource(W, Point(0, -0.5), 1)
delta.apply(b)
delta = PointSource(W, Point(0, 0.5), 1)
delta.apply(b)

w = Function(W)
solve(A, w.vector(), b)
(u, c) = w.split()
vtkfile = File('neumannDirac/numerical.pvd')
vtkfile << u
