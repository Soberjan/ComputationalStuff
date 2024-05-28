import matplotlib.pyplot as plt
import numpy as np
from fenics import *

T = 10
num_steps = 100
dt = T / num_steps

Du = 5.0
L = 15

#Меш
mesh = IntervalMesh(100, 0.0, L)
OX = mesh.coordinates()
V = FunctionSpace(mesh, "Lagrange", 1)
figure, axis = plt.subplots(2, 1)
axis[0].set_title("Analytical")
axis[1].set_title("Numerical")

#Аналитическое решение
def U(x, t):
    return np.cos(np.pi*x/L) * 5*np.exp(-Du*np.pi**2*t/L**2)

t = 0
dt2 = 5
while t <= 10:
    uAnal = U(OX, t) #AAAARGH
    axis[0].plot(OX, uAnal, label="t="+str(t))
    t += dt2

#Численное решение
u_n = Function(V)
u_0 = Expression("5*cos(3.14*x[0]/L)", degree=0, L=L)
#u_0 = Expression("2", degree=0, L=L)
u_n = interpolate(u_0, V)
u_n_vec = u_n.compute_vertex_values(mesh)
axis[1].plot(OX, u_n_vec, label="t="+str(0.0))

u = Function(V)
v = TestFunction(V)

Du = Constant(Du)
k = Constant(dt)
F = (u-u_n)/dt*v*dx + Du*dot(grad(u), grad(v))*dx
#F = (u-u_n)/k*v*dx + Du*u.dx(0)*v.dx(0)*dx

t = 0
n = 0
for n in range(num_steps):
    t += dt
    n += 1
    solve(F == 0, u)
    u_n.assign(u)
    if n % 50 == 0:
        u_n_vec = u_n.compute_vertex_values(mesh)
        axis[1].plot(OX, u_n_vec, label="t="+str(round(t)))

axis[0].legend()
axis[1].legend()
plt.show()