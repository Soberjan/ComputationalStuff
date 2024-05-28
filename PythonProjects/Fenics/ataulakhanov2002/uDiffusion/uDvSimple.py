import matplotlib.pyplot as plt
from fenics import *
import numpy as np

T = 10
num_steps = 100
dt = T / num_steps

Du = 5.0
a = 1.0
eps = 0.01

L = 15.0
mesh = IntervalMesh(num_steps, 0.0, L)
OX = mesh.coordinates()
P1 = FiniteElement("Lagrange", interval, 1)
element = MixedElement([P1, P1])
V = FunctionSpace(mesh, element)

v_1, v_2 = TestFunctions(V)

u = Function(V)
u_1, u_2 = split(u)
u_0 = Expression(("5*cos(3.14*x[0]/L)", "0"), degree=0, L=L)
u_n = Function(V)
u_n = interpolate(u_0, V)
u_n_vec = u_n.compute_vertex_values(mesh)
u_n1, u_n2 = split(u_n)
u_n1_vec, u_n2_vec = np.split(u_n_vec, 2)
figure, axis = plt.subplots(2, 1)
axis[0].plot(OX, u_n1_vec, label=0.0)
axis[1].plot(OX, u_n2_vec, label=0.0)

Du = Constant(Du)
a = Constant(a)
k = Constant(dt)
eps = Constant(eps)

F = (u_1-u_n1)/dt*v_1*dx + Du*dot(grad(u_1), grad(v_1))*dx \
    + (u_2-u_n2)/dt*v_2*dx - eps*a*u_1*v_2*dx + eps*a*u_2*v_2*dx

t = 0
n = 0
for n in range(num_steps):
    t += dt
    n+=1
    solve(F == 0, u)
    _u_1, _u_2 = u.split()
    u_vec = u.compute_vertex_values(mesh)
    _u_1_vec, _u_2_vec = np.split(u_vec, 2)
    if n % 50 == 0:
        axis[0].plot(OX, _u_1_vec, label=t)
        axis[1].plot(OX, _u_2_vec, label=t)
    u_n.assign(u)

plt.legend()
plt.show()