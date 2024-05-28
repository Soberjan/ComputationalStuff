import matplotlib.pyplot as plt
from fenics import *

T = 10
num_steps = 100
dt = T / num_steps

Du = 5
L = 15
n = 0.4

mesh = IntervalMesh(100, 0.0, L)
V = FunctionSpace(mesh, "Lagrange", 1)

u_n = Function(V)
u_0 = Expression("cos(3.14*x[0]/L)", degree=0, L=L)
#u_0 = Expression("x[0]", degree=0, L=L)
u_n = interpolate(u_0, V)
plt.figure()
#plt.ylim(-200, 200)
plot(u_n, label=0.0)

u = Function(V)
v = TestFunction(V)

Du = Constant(Du)
k = Constant(dt)
F = (u-u_n)/dt*v*dx + Du*dot(grad(u), grad(v))*dx + u*(u-1)*(u-n)*v*dx
#F = (u-u_n)/k*v*dx + Du*u.dx(0)*v.dx(0)*dx

t = 0
num = 0
for n in range(num_steps):
    t += dt
    num += 1
    solve(F == 0, u)
    u_n.assign(u)
    if num % 10 == 0:
        plot(u_n, label=t)

plt.legend()
plt.show()
