from fenics import *
import numpy as np
import matplotlib.pyplot as plt


def solve_dirichlet():
    # print(n)
    # mesh = UnitSquareMesh(5, 5)
    mesh = Mesh("/home/soberjan/PythonProjects/Fenics/EEG/Debugging/Meshes/donut.xml")
    
    V = FunctionSpace(mesh, 'P', 1)

    u_D = Expression('1 + x[0]*x[0] + 2*x[1]*x[1]', degree=2)

    def boundary(x, on_boundary):
        return on_boundary

    bc = DirichletBC(V, u_D, boundary)

    u = TrialFunction(V)
    v = TestFunction(V)
    f = Constant(-6.0)
    a = dot(grad(u), grad(v)) * dx
    L = f * v * dx

    u = Function(V)
    solve(a == L, u, bc)

    vtkfile = File('dirichlet/numerical.pvd')
    vtkfile << u
    
    # analytical = Function(V)
    # analytical.interpolate(u_D)
    # vtkfile1 = File('dirichlet/analytical.pvd')
    # vtkfile1 << analytical



    # Compute error in L2 norm
    error_L2 = errornorm(u_D, u, 'L2')
    # Compute maximum error at vertices
    vertex_values_u_D = u_D.compute_vertex_values(mesh)
    vertex_values_u = u.compute_vertex_values(mesh)
    error_max = np.max(np.abs(vertex_values_u_D - vertex_values_u))
    # Print errors
    #print('error_L2 =', error_L2)
    #print('error_max =', error_max)
    return error_L2

print(solve_dirichlet())

# OX = np.array([2, 4, 8, 16, 32, 64, 128, 256])
# OY = []
# for n in OX:
#     OY.append(solve_dirichlet(n))

# plt.plot(OX**2, OY)

# plt.xlabel('number of nodes')
# plt.ylabel('L2 norm')
# plt.show()
