from dolfin import *
import mpi4py.MPI

k1 = 9
k2 = 10

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
f1 = Expression('x[0] < -0.5 + eps && x[0] > -0.5 - eps && x[1] < -0.5 + eps && x[1] > -0.5 - eps ? sin(k1 * t) : 0', k1=k1, eps=eps, t=t, degree=1)
f2 = Expression('x[0] < 0 + eps && x[0] > 0 - eps && x[1] < 0.5 + eps && x[1] > 0.5 - eps ? sin(k2 * t) : 0', k2=k2, eps=eps, t=t, degree=1)
f3 = Expression('x[0] < 0.5 + eps && x[0] > 0.5 - eps && x[1] < -0.5 + eps && x[1] > -0.5 - eps ? sin(k1 * t) + sin(k2 * t) : 0', k1=k1, k2=k2, eps=eps, t=t, degree=1)

g = Expression("0", degree=2)
w = Function(W)
a = (inner(grad(u), grad(v)) + c * v + u * d) * dx
L = (f1 + f2 - f3) * v * dx

file = File("RotatingRegime/RotatingRegime.pvd")
dt = 0.01
while t <= 2:
    f1.t = t
    f2.t = t
    f3.t = t
    solve(a == L, w)
    (u, c) = w.split()
    file << (u, t)
    print(t)
    t += dt

