from fenics import *
import matplotlib.pyplot as plt

T = 2.0
num_steps = 100
dt = T / num_steps
eps = 0.01
K = 10.0

mesh = UnitSquareMesh(10, 10)
P1 = FiniteElement('Lagrange', triangle, 1)
element = MixedElement([P1, P1, P1])
V = FunctionSpace(mesh, element)

v_1, v_2, v_3 = TestFunctions(V)

u = Function(V)
u_n = Function(V)
u_1, u_2, u_3 = split(u)
u_n1, u_n2, u_n3 = split(u_n)

f_1 = Expression("pow(x[0]-0.1,2)+pow(x[1]-0.1,2)<0.05*0.05 ? 0.1 : 0", degree=1)
f_2 = Expression("pow(x[0]-0.1,2)+pow(x[1]-0.3,2)<0.05*0.05 ? 0.1 : 0", degree=1)
f_3 = Expression("pow(x[0]-0.1,2)+pow(x[1]-0.3,2)<0.05*0.05 ? 0.1 : 0", degree=1)

k = Constant(dt)
K = Constant(K)
eps = Constant(eps)

F = ((u_1 - u_n1) / k)*v_1*dx + grad(u_1)*v_1*dx\
    + eps*dot(grad(u_1), grad(v_1))*dx + K*u_1*u_2*v_1*dx \
    + ((u_2 - u_n2) / k)*v_2*dx + dot(grad(u_2), grad(u_2))*v_2*dx \
    + eps*dot(grad(u_2), grad(v_2))*dx + K*u_1*u_2*v_2*dx \
    + ((u_3 - u_n3) / k)*v_3*dx + dot(grad(u_3), grad(u_3))*v_3*dx \
    + eps*dot(grad(u_3), grad(v_3))*dx - K*u_1*u_2*v_3*dx + K*u_3*v_3*dx \
    - f_1*v_1*dx - f_2*v_2*dx - f_3*v_3*dx

vtkfile_u_1 = File("reaction_system/u_1.pvd")
vtkfile_u_2 = File("reaction_system/u_2.pvd")
vtkfile_u_3 = File("reaction_system/u_3.pvd")

t = 0
for n in range(num_steps):
    t += dt

    solve(F == 0, u)

    _u_1, _u_2, _u_3 = u.split()
    vtkfile_u_1 << (_u_1, t)
    vtkfile_u_2 << (_u_2, t)
    vtkfile_u_3 << (_u_3, t)

    u_n.assign(u)
