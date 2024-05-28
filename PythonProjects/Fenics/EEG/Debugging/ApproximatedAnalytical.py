from dolfin import *
import numpy as np
import mpi4py.MPI
import matplotlib.pyplot as plt
import meshio


def solve_dirac(n, save_res=False):
    comm = mpi4py.MPI.COMM_WORLD
    mesh = UnitDiscMesh.create(comm, n, 2, 2)

    # Build function space with Lagrange multiplier
    P1 = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
    R = FiniteElement("Real", mesh.ufl_cell(), 0)
    W = FunctionSpace(mesh, P1 * R)
    V = FunctionSpace(mesh, P1)

    eps = 1/n
    point1 = [0, 0]
    f1 = Expression('x[0] < x_i+eps && x[0] > x_i-eps && x[1] < y_i+eps && x[1] > y_i-eps ? 1 / (4 * eps * eps) : 0',
                    eps=eps, degree=1, x_i=point1[0], y_i=point1[1])
    # point2 = [0.0, -0.5]
    # f2 = Expression('x[0] < x_i+eps && x[0] > x_i-eps && x[1] < y_i+eps && x[1] > y_i-eps ? 1 / (4 * eps * eps) : 0',
    #                 eps=eps, degree=1, x_i=point2[0], y_i=point2[1])
    u_anal = Expression('x[0] == x_1 && x[1] == y_1 ? 100 : -1/(4*pi*pow(pow(x[0]-x_1, 2)+pow(x[1]-y_1, 2), 0.5))',
                        degree=1, x_1=point1[0], y_1=point1[1], pi=np.pi)

    # Define variational problem
    (u, c) = TrialFunctions(W)
    (v, d) = TestFunctions(W)

    a = (inner(grad(u), grad(v)) + c * v + u * d) * dx
    L = (f1) * v * dx

    w = Function(W)
    solve(a == L, w)
    (u, c) = w.split()

    if save_res:
        vtkfile = File('approximatedAnalytical/numerical.pvd')
        vtkfile << u

        k = Function(V)
        k.interpolate(u_anal)
        vtkfile = File('approximatedAnalytical/analytical.pvd')
        vtkfile << k

    error_L2 = errornorm(u_anal, u, 'l2')
    return error_L2

# solve_dirac(64, True)

OX = np.array([4, 8, 16, 32, 64, 128])
OY = []
for n in OX:
    OY.append(solve_dirac(n))
plt.plot(OX**2, OY)
plt.xlabel('number of nodes')
plt.ylabel('L2 error')
plt.show()

#mesh = Mesh("EEG/Debugging/Meshes/donut.xml")
