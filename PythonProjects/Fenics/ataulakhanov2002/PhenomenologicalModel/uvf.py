import matplotlib.pyplot as plt
import numpy as np
from fenics import *

T = 24
num_steps = 200
dt = T / num_steps

k1 = 2.0
k2 = 100
k3 = 0.05
k4 = 5
k5 = 0.0015
k6 = 5
k7 = 0.05
k8 = 0.35
k9 = 2.8
Du = 0.001
Dv = 0.005

L = 1.5
mesh = IntervalMesh(1000, 0.0, L)
OX = mesh.coordinates()
P1 = FiniteElement("Lagrange", interval, 1)
element = MixedElement([P1, P1, P1])
V = FunctionSpace(mesh, element)

v_1, v_2, v_3 = TestFunctions(V)

u = Function(V)
u_1, u_2, u_3 = split(u)
#u_0 = Expression(("0.3*cos(3.14*x[0]/L)", "0.07"), degree=0, L=L)
tol = 1E-14
#u_0 = Expression(("(x[0] >= -tol) && (x[0] <= 0.1 + tol) ? 0.2 : 0.09395", "-0.027", "0"), degree=0, L=L, tol=tol)
u_0 = Expression(("0", "10", "0"), degree=0, L=L, tol=tol)
u_n = Function(V)
u_n = interpolate(u_0, V)
u_n_vec = u_n.compute_vertex_values(mesh)
u_n1, u_n2, u_n3 = split(u_n)
u_n1_vec, u_n2_vec, u_n3_vec = np.split(u_n_vec, 3)
figure, axis = plt.subplots(3, 1)
axis[0].plot(OX, u_n1_vec, label=0.0)
axis[1].plot(OX, u_n2_vec, label=0.0)
axis[2].plot(OX, u_n3_vec, label=0.0)

k1 = Constant(k1)
k2 = Constant(k2)
k3 = Constant(k3)
k4 = Constant(k4)
k5 = Constant(k5)
k6 = Constant(k6)
k7 = Constant(k7)
k8 = Constant(k8)
k9 = Constant(k9)
Du = Constant(Du)
Dv = Constant(Dv)
dt = Constant(dt)

F = (u_1-u_n1)/dt*v_1*dx - k1*u_1**2/(u_1+k2)*v_1*dx + k3*u_1*v_1*dx + k4*u_1*u_2*v_1*dx + Du*dot(grad(u_1), grad(v_1))*dx \
    + (u_2-u_n2)/dt*v_2*dx - k5*u_1*(1-u_2/k6)*(1+(u_2/k7)**2)*v_2*dx + k8*u_2*v_2*dx + Dv*dot(grad(u_2), grad(v_2))*dx \
    + (u_3-u_n3)/dt*v_3*dx - k9*u_1*v_3*dx

# F = (u_1-u_n1)/dt*v_1*dx - k1*(u_1**2)/(u_1+k2)*v_1*dx + k3*u_1*v_1*dx + k4*u_1*u_2*v_1*dx + Du*dot(grad(u_1), grad(v_1))*dx \
#     + (u_2-u_n2)/dt*v_2*dx - k5*u_1*(1-u_2/k6)*(1+(u_2/k7)**2)*v_2*dx + k8*u_2*v_2*dx + Dv*dot(grad(u_2), grad(v_2))*dx \
#     + (u_3-u_n3)/dt*v_3*dx - k9*u_1*v_3*dx

#F = (u_1-u_n1)/dt*v_1*dx - u_1/(u_1+k2)*v_1*dx + Du*dot(grad(u_1), grad(v_1))*dx


t = 0
n = 0
for n in range(num_steps):
    t += dt
    n+=1
    solve(F == 0, u)
    _u_1, _u_2, _u_3 = u.split()
    u_vec = u.compute_vertex_values(mesh)
    _u_1_vec, _u_2_vec, _u_3_vec = np.split(u_vec, 3)
    if n % 100 == 0:
        axis[0].plot(OX, _u_1_vec, label=t)
        axis[1].plot(OX, _u_2_vec, label=t)
        axis[2].plot(OX, _u_3_vec, label=t)
    u_n.assign(u)
    print(n)
# #plt.legend()
plt.show()