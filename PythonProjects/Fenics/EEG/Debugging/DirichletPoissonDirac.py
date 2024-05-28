from fenics import *
import numpy as np

mesh = UnitSquareMesh(20, 20)
V = FunctionSpace(mesh, 'P', 1)

pi = np.pi
u_D = Expression('x[0] == 0 && x[1] == 0 ? 0 : '
                 '1 / (2 * pi) * std::log(sqrt(x[0] * x[0] + x[1] * x[1]))', degree=1, pi=pi)

def boundary(x, on_boundary):
    return on_boundary
bc = DirichletBC(V, u_D, boundary)

u = TrialFunction(V)
v = TestFunction(V)

a = dot(grad(u), grad(v)) * dx
L = Constant(0) * v * dx
A, b = assemble_system(a, L, bc)

delta = PointSource(V, Point(0.5, 0.5), 1)
delta.apply(b)

u = Function(V)
solve(A, u.vector(), b)

vtkfile = File('dirichletDirac/numerical.pvd')
vtkfile << u

analytical = Function(V)
analytical.interpolate(u_D)
vtkfile1 = File('dirichletDirac/analytical.pvd')
vtkfile1 << analytical

# Compute error in L2 norm
error_L2 = errornorm(u_D, u, 'L2')
print('error_L2 =', error_L2)
