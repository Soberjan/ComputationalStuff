from dolfin import *
import mpi4py.MPI
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection


def plot_mesh(coordinates, n, m, type):
    coordinates = np.transpose(coordinates)
    if type == 'Circle':
        plt.scatter(coordinates[0], coordinates[1])
    if type == 'Square':
        x = np.reshape(coordinates[0], (n+1, m+1))
        y = np.reshape(coordinates[1], (n+1, m+1))
        plt.scatter(coordinates[0], coordinates[1])
        segs1 = np.stack((x, y), axis=2)
        segs2 = segs1.transpose(1, 0, 2)
        plt.gca().add_collection(LineCollection(segs1))
        plt.gca().add_collection(LineCollection(segs2))


# comm = mpi4py.MPI.COMM_WORLD
# mesh = UnitDiscMesh.create(comm, 1, 1, 2)
n = 10
m = 10
mesh = UnitSquareMesh(n, m)
V = FunctionSpace(mesh, "Lagrange", 1)

plot_mesh(mesh.coordinates(), n, m, 'Square')
#plt.show()

def boundary(x):
    return x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS
u0 = Constant(0.0)
bc = DirichletBC(V, u0, boundary)

u = TrialFunction(V)
v = TestFunction(V)
#f = Expression("10*exp(-(pow(x[0] - 0.5, 2) + pow(x[1] - 0.5, 2)) / 0.02)", degree=2)
g = Expression("sin(5*x[0])", degree=1)
a = inner(grad(u), grad(v))*dx
L = g*v*dx

u = Function(V)
solve(a == L, u)

vertex_values = u.compute_vertex_values()
print(vertex_values)
