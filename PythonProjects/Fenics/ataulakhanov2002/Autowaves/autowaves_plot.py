import matplotlib.pyplot as plt
from fenics import *
import numpy as np

T = 1500
num_steps = 1000
dt = T / num_steps

Du = 0.001
eps = 0.001
a = 4
b = -0.55
n = 0.4
L = 15
mesh = IntervalMesh(450, 0.0, L)
OX = mesh.coordinates()
P1 = FiniteElement("Lagrange", interval, 1)
element = MixedElement([P1, P1])
V = FunctionSpace(mesh, element)

v_1, v_2 = TestFunctions(V)

u = Function(V)
u_1, u_2 = split(u)
tol = 1E-14
u_0 = Expression(("0.804", "(x[0] <= 0.1) ? 0.08 : 0.064"), degree=0, L=L, tol=tol)
u_n = Function(V)
u_n = interpolate(u_0, V)
u_n1, u_n2 = split(u_n)

u_n_vec = u_n.compute_vertex_values(mesh)
u_n1_vec, u_n2_vec = np.split(u_n_vec, 2)
figure, axis = plt.subplots(2, 1)
axis[0].plot(OX, u_n1_vec, label=0.0)
axis[1].plot(OX, u_n2_vec, label=0.0)
axis[0].set_ylim(-0.25, 1.25)
axis[1].set_ylim(-0.1, 0.2)

Du = Constant(Du)
a = Constant(a)
k = Constant(dt)
eps = Constant(eps)
b = Constant(b)
n = Constant(n)

F = (u_1-u_n1)/dt*v_1*dx + (u_1*(u_1-1)*(u_1-n)+u_2)*v_1*dx + Du*dot(grad(u_1), grad(v_1))*dx +\
    (u_2-u_n2)/dt*v_2*dx - eps*(u_1-a*u_2+b)*v_2*dx

t = 0
n = 0

for n in range(num_steps):
    print(n)
    t += dt
    solve(F == 0, u)
    _u_1, _u_2 = u.split()
    if n % 200 == 0:
        u_vec = u.compute_vertex_values(mesh)
        _u_1_vec, _u_2_vec = np.split(u_vec, 2)
        axis[0].plot(OX, _u_1_vec, label = t)
        axis[1].plot(OX, _u_2_vec, label = t)
    n += 1
    u_n.assign(u)

plt.legend()
plt.show()