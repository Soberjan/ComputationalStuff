from dolfin import *
import numpy as np

mesh = UnitSquareMesh(16,8)
V = FunctionSpace(mesh,'CG',1)

C1 = Function(V)

class InitialCondition(Expression):
    def eval_cell(self, value, x, ufc_cell):
        if x[0] <= 0.5:
            value[0] = 1.0
        else:
            value[0] = 0.0

C1.interpolate(InitialCondition())
plot(C1, interactive=True)