from fenics import *

T = 5.0
num_steps = 100
dt = T / num_steps
a = 1
b = 1
tau = 1

mesh = UnitSquareMesh(100, 100)
P1 = FiniteElement('Lagrange', triangle, 1)
element = MixedElement([P1, P1])
V = FunctionSpace(mesh, element)

v_1, v_2 = TestFunctions(V)

u = Function(V)
u_1, u_2 = split(u)

u_0 = Expression(("sin(x[0])", "cos(x[0]*x[1])"), degree=1)
u_n = Function(V)
u_n = project(u_0, V)
u_n1, u_n2 = split(u_n)

# Как распространить граничные условия на все функции?
def boundary(x, on_boundary):
    return on_boundary
u_D = Expression("2 * sin(x[0])", degree=1)
bc = DirichletBC(V.sub(1), u_D, boundary)

k = Constant(dt)
a = Constant(a)
b = Constant(b)
tau = Constant(tau)

F = ((u_1-u_n1)/k)*v_1*dx - u_1*v_1*dx + ((u_1**3)/3.0)*v_1*dx + u_2*v_1*dx \
    + ((u_2-u_n2)/k*tau)*v_2*dx - (u_1+a)*v_2*dx + (b*u_2)*v_2*dx

vtkfile_u_1 = File("fitzhugh/u_1.pvd")
vtkfile_u_2 = File("fitzhugh/u_2.pvd")

t = 0
for n in range(num_steps):
    t += dt

    solve(F == 0, u, bc)

    _u_1, _u_2= u.split()
    vtkfile_u_1 << (_u_1, t)
    vtkfile_u_2 << (_u_2, t)

    u_n.assign(u)
