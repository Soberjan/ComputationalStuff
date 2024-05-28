from fenics import *
import numpy as np
import matplotlib.pyplot as plt

def solve_dirac(n):
    mesh = UnitSquareMesh(n, n)
    V = FunctionSpace(mesh, 'P', 1)

    pi = np.pi
    u_D = Expression('x[0] == 0 && x[1] == 0 ? 0 : '
                     '1 / (2 * pi) * std::log(sqrt(x[0] * x[0] + x[1] * x[1]))', degree=1, pi=pi)

    def boundary(x, on_boundary):
        return on_boundary
    bc = DirichletBC(V, u_D, boundary)

    u = TrialFunction(V)
    v = TestFunction(V)

    eps = 1
    f1 = Expression('x[0] < eps && x[0] > eps && x[1] < eps && x[1] > - eps ? 1 / (4 * eps * eps) : 0', eps=eps, degree=1)

    a = dot(grad(u), grad(v)) * dx
    L = f1 * v * dx

    u = Function(V)
    solve(a == L, u, bc)

    # vtkfile = File('dirichletDirac/numerical_expr.pvd')
    # vtkfile << u

    analytical = Function(V)
    analytical.interpolate(u_D)
    # vtkfile1 = File('dirichletDirac/analytical.pvd')
    # vtkfile1 << analytical

    # Compute error in L2 norm
    # u_vec = u.compute_vertex_values()
    # anal_vec = analytical.compute_vertex_values()
    # u_dif = u_vec - anal_vec
    # return (sum(u_dif**2)**0.5)
    error_L2 = errornorm(u_D, u, 'L2')
    return error_L2

OX = np.array([2, 4, 8, 16, 32, 64, 128, 256])
OY = []
for n in OX:
    OY.append(solve_dirac(n))

plt.plot(OX**2, OY)
plt.xlabel('number of nodes')
plt.ylabel('L2 error')
plt.show()
