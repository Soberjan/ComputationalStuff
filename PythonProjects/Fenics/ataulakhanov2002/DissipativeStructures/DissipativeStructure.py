import matplotlib.pyplot as plt
from fenics import *
#на границе производная u,v = 0
#попробовать смоделировать с Du=1, Dv=0
T = 1000
num_steps = 1000
dt = T / num_steps

Du = 0.001
Dv = 0.005
n = 0.4
a = 9.3995
b = -0.2637
eps = 0.02

L = 3.0
mesh = IntervalMesh(30, 0.0, L)
P1 = FiniteElement("Lagrange", interval, 1)
element = MixedElement([P1, P1])
V = FunctionSpace(mesh, element)

v_1, v_2 = TestFunctions(V)

u = Function(V)
u_1, u_2 = split(u)
tol = 1E-14
u_0 = Expression(("(x[0] >= (L/2 - 0.1 - tol) && x[0] <= (L/2 + 0.1 + tol)) ? 0.742 + sin(x[0]) : 0.742", "0.065"), degree=0, tol=tol, L=L)
u_n = Function(V)
u_n = interpolate(u_0, V)
u_n1, u_n2 = split(u_n)
plt.figure()
plot(u_n1, label="t=0.0")

Du = Constant(Du)
Dv = Constant(Dv)
n = Constant(n)
a = Constant(a)
b = Constant(b)
k = Constant(dt)
eps = Constant(eps)

F = ((u_1-u_n1)/k)*v_1*dx + Du*u_1.dx(0)*v_1.dx(0)*dx + u_1*(u_1-1)*(u_1-n)*v_1*dx + u_2*v_1*dx \
    + ((u_2-u_n2)/k)*v_2*dx - eps*(u_1-a*u_2+b)*v_2*dx + Dv*u_2.dx(0)*v_2.dx(0)*dx
# F = (u_1-u_n1)/k*v_1*dx + u_1*(u_1-1)*(u_1-n)*v_1*dx + u_2*v_1*dx \
#     + (u_2-u_n2)/k*v_2*dx - eps * (u_1-a*u_2+b)*v_2*dx

# vtkfile_u_2 = File("diffusion_u/v.pvd")
# vtkfile_u_1 = File("diffusion_u/u.pvd")

t = 0
for n in range(num_steps):
    t += dt

    solve(F == 0, u)

    _u_1, _u_2 = u.split()
    # vtkfile_u_1 << (_u_1, t)
    # vtkfile_u_2 << (_u_2, t)
    uselessVariable = int(t) % 100
    print(t)
    if uselessVariable == 0:
        plot(_u_1, label=t)
    u_n.assign(u)

plt.legend()
plt.show()