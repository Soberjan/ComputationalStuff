from fenics import *

T = 50
num_steps = 1000
dt = T / num_steps
a = 0.03125
c_1 = 9
c_2 = 0.02
c_3 = 5
eps = 0.1
D1 = D2 = 1

L = 400.0
H = 32.0
mesh = RectangleMesh(Point(0.0, 0.0), Point(L, H), 100, 32, diagonal="right")
P1 = FiniteElement("Lagrange", triangle, 1)
element = MixedElement([P1, P1])
V = FunctionSpace(mesh, element)

v_1, v_2 = TestFunctions(V)

u = Function(V)
u_1, u_2 = split(u)

tol = 1E-14
u_0 = Expression(("(x[0] >= (L/2 - 10 - tol) && x[0] <= (L/2 + 10 + tol)) ? 1 : 0", "0"), degree=0, tol=tol, L=L)
u_n = Function(V)
u_n = interpolate(u_0, V)
u_n1, u_n2 = split(u_n)

# def boundary(x, on_boundary):
#     return on_boundary
# u_D_1 = Expression("0", degree=1)
# bc = DirichletBC(V.sub(0), u_D_1, boundary)

k = Constant(dt)
c_1 = Constant(c_1)
c_2 = Constant(c_2)
c_3 = Constant(c_3)
eps = Constant(eps)
D1 = Constant(D1)
D2 = Constant(D2)

#g = Expression("sin(x[1])*3.14*2", degree=1)
V_y = Expression("a*x[1]*(H-x[1])", degree=2, a=a, H=H)
g = Expression("0", degree=0)

F = ((u_1-u_n1)/k)*v_1*dx - c_1*u_1*(u_1-c_2)*(1-u_1)*v_1*dx + u_2*v_1*dx + D1*dot(grad(u_1), grad(v_1))*dx + D1*g*v_1*ds + V_y*u_1.dx(0)*v_1*dx\
    + ((u_2-u_n2)/k)*v_2*dx - eps*(c_3*u_1-u_2)*v_2*dx + D2*dot(grad(u_2), grad(v_2))*dx + D2*g*v_2*ds + V_y*u_2.dx(0)*v_2*dx

vtkfile_u_1 = File("fh-flow_T2/u_1.pvd")
vtkfile_u_2 = File("fh-flow_T2/u_2.pvd")

t = 0
for n in range(num_steps):
    t += dt

    solve(F == 0, u)

    _u_1, _u_2 = u.split()
    vtkfile_u_1 << (_u_1, t)
    vtkfile_u_2 << (_u_2, t)

    u_n.assign(u)
