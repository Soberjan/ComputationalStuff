from fenics import *
import matplotlib.pyplot as plt
import numpy as np

xmax = np.pi
ymax = np.pi
mesh = RectangleMesh(Point(0.0, 0.0), Point(xmax, ymax), 100, 100, diagonal="right")
V = FunctionSpace(mesh, 'P', 1)

u_Left = Expression("0", degree=0)
def boundary_Left(x, on_boundary):
    tol = 1E-14
    return on_boundary and near(x[0], 0, tol)
bc_Left = DirichletBC(V, u_Left, boundary_Left)

u_Right = Expression("0", degree=0)
def boundary_Right(x, on_boundary):
    tol = 1E-14
    return on_boundary and near(x[0], xmax, tol)
bc_Right = DirichletBC(V, u_Right, boundary_Right)

bcs = [bc_Left, bc_Right]

u = TrialFunction(V)
v = TestFunction(V)
f = Expression("0", degree=1)
a = dot(grad(u), grad(v))*dx
tol = 1E-14
#g = Expression("x[1] <= 3.14 + tol ? sin(3*x[0]) : 0", degree=1, tol=tol)
g = Expression("x[1]/3.14*sin(3*x[0])", degree=1)
L = f*v*dx - g*v*ds
#u_0 = Expression(("(x[0] >= (L/2 - 10 - tol) && x[0] <= (L/2 + 10 + tol)) ? 1 : 0", "0"), degree=0, tol=tol, L=L)


vtkfile_u = File("Task/u.pvd")
u = Function(V)
solve(a == L, u, bcs)
plot(u)
plot(mesh)
vtkfile_u << (u, 0)

plt.show()