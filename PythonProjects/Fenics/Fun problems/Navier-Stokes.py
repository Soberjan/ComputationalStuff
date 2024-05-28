from fenics import *
from mshr import *
import numpy as np

T = 5.0
num_steps = 5000
dt = T / num_steps
mu = 0.001
rho = 1

channel = Rectangle(Point(0, 0), Point(2.2, 0.41))
cylinder = Circle(Point(0.2, 0.2), 0.05)
domain = channel - cylinder
mesh = generate_mesh(domain, 64)

V = VectorFunctionSpace(mesh, ’P’, 2)
Q = FunctionSpace(mesh, ’P’, 1)