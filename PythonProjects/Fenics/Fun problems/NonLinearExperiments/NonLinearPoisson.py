from fenics import *
import matplotlib.pyplot as plt

mesh = UnitSquareMesh(10, 10)
V = FunctionSpace(mesh, "Lagrange", 1)

tol = 1E-14
def left_boundary(x, on_boundary):
    return on_boundary and abs(x[0]) < tol

def right_boundary(x, on_boundary):
    return on_boundary and abs(x[0]-1) < tol

Gamma_0 = DirichletBC(V, Constant(0.0), left_boundary)
Gamma_1 = DirichletBC(V, Constant(1.0), right_boundary)
bcs = [Gamma_0, Gamma_1]
m = 2
def q(u):
    return (1+u)**m

def Dq(u):
    return m*(1+u)**(m-1)

u = TrialFunction(V)
v = TestFunction(V)
u_ = Function(V)      # the most recently computed solution
F = inner(q(u)*nabla_grad(u), nabla_grad(v))*dx
F = action(F, u_)

J = derivative(F, u_, u)   # Gateaux derivative in dir. of u

problem = NonlinearVariationalProblem(F, u_, bcs, J)
solver = NonlinearVariationalSolver(problem)

prm = solver.parameters
prm['newton_solver']['absolute_tolerance'] = 1E-8
prm['newton_solver']['relative_tolerance'] = 1E-7
prm['newton_solver']['maximum_iterations'] = 25
prm['newton_solver']['relaxation_parameter'] = 1.0
# if iterative_solver:
#     prm['linear_solver'] = 'gmres'
#     prm['preconditioner'] = 'ilu'
#     prm['krylov_solver']['absolute_tolerance'] = 1E-9
#     prm['krylov_solver']['relative_tolerance'] = 1E-7
#     prm['krylov_solver']['maximum_iterations'] = 1000
#     prm['krylov_solver']['gmres']['restart'] = 40
#     prm['krylov_solver']['preconditioner']['ilu']['fill_level'] = 0
#set_log_level(PROGRESS)

solver.solve()

plot(u_)
#plot(u)
plt.show()