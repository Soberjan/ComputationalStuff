import matplotlib.pyplot as plt
import numpy as np
from fenics import *

T = 1000
num_steps = 5000
dt = T / num_steps

Du = 0.001
Dv = 0.07
a = 4
b = -0.48
eps = 0.01
n = 0.4

# Du = 0.001
# Dv = 0.005
# a = 9.3995
# b = -0.405
# eps = 0.02
# n = 0.4

L = 3.0
mesh = IntervalMesh(1000, 0.0, L)
OX = mesh.coordinates()
P1 = FiniteElement("Lagrange", interval, 1)
element = MixedElement([P1, P1])
V = FunctionSpace(mesh, element)

v_1, v_2 = TestFunctions(V)

u = Function(V)
u_1, u_2 = split(u)
#u_0 = Expression(("0.3*cos(3.14*x[0]/L)", "0.07"), degree=0, L=L)
tol = 1E-14
u_0 = Expression(("(x[0] >= L/2-0.1-tol) && (x[0] <= L/2+0.1+tol) ? 0.09396 : 0.09395", "-7.4"), degree=0, L=L, tol=tol)
u_n = Function(V)
u_n = interpolate(u_0, V)
u_n_vec = u_n.compute_vertex_values(mesh)
u_n1, u_n2 = split(u_n)
u_n1_vec, u_n2_vec = np.split(u_n_vec, 2)
figure, axis = plt.subplots(2, 1)
axis[0].plot(OX, u_n1_vec, label=0.0)
axis[0].set_title('U')
axis[1].plot(OX, u_n2_vec, label=0.0)
axis[1].set_title('V')

Du = Constant(Du)
Dv = Constant(Dv)
a = Constant(a)
k = Constant(dt)
eps = Constant(eps)
b = Constant(b)
n = Constant(n)

F = (u_1-u_n1)/dt*v_1*dx + Du*dot(grad(u_1), grad(v_1))*dx\
    + (u_2-u_n2)/dt*v_2*dx - eps*a*u_1*v_2*dx + eps*a*u_2*v_2*dx - eps*b*v_2*dx + Dv*dot(grad(u_2), grad(v_2))*dx

# F = (u_1-u_n1)/dt*v_1*dx + Du*dot(grad(u_1), grad(v_1))*dx + u_1**3*v_1*dx \
#     + (u_2-u_n2)/dt*v_2*dx - eps*a*u_1*v_2*dx + eps*a*u_2*v_2*dx - eps*b*v_2*dx + Dv*dot(grad(u_2), grad(v_2))*dx

t = 0
n = 0
for n in range(num_steps):
    t += dt
    n+=1
    solve(F == 0, u)
    _u_1, _u_2 = u.split()
    u_vec = u.compute_vertex_values(mesh)
    _u_1_vec, _u_2_vec = np.split(u_vec, 2)
    if n % 100 == 0:
        axis[0].plot(OX, _u_1_vec, label=t)
        axis[1].plot(OX, _u_2_vec, label=t)
    u_n.assign(u)
    print(n)
#plt.legend()
plt.show()