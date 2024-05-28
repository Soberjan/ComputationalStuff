from dolfin import *
import mpi4py.MPI

k1 = 1
k2 = 1

comm = mpi4py.MPI.COMM_WORLD
mesh = UnitDiscMesh.create(comm, 25, 2, 2)

# Build function space with Lagrange multiplier
P1 = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
R = FiniteElement("Real", mesh.ufl_cell(), 0)
W = FunctionSpace(mesh, P1 * R)

# Define variational problem
(u, c) = TrialFunction(W)
(v, d) = TestFunctions(W)

a = (inner(grad(u), grad(v)) + c * v + u * d) * dx
L = Constant(0) * v * dx
A, b = assemble_system(a, L)

file = File("StandingWaves/StandingWaves3Sources22.pvd")

w = Function(W)
t = 0
dt = 0.01
while t <= 5:
    delta = PointSource(W, Point(-0.5, -0.5), sin(k1*t))
    delta.apply(b)
    delta = PointSource(W, Point(0, 0.5), sin(k2*t))
    delta.apply(b)
    delta = PointSource(W, Point(0.5, -0.5), -(sin(k1*t)+sin(k2*t)))
    delta.apply(b)

    solve(A, w.vector(), b)
    (u, c) = w.split()
    file << (u, t)
    print(t)
    t += dt

