import numpy as np
import matplotlib.pyplot as plt

def u(x, y):
    return -(np.sin(x)*(-0.4)*np.cos(2*y) + np.sin(3 * x) * ((np.exp(3 * y) + np.exp(-3 * y)) / (3 * (np.exp(3 * np.pi) - np.exp(-3 * np.pi)))))

x = np.linspace(0, np.pi, 100)
y = np.linspace(0, np.pi, 100)

X, Y = np.meshgrid(x, y)
nu = u(X, Y)
plt.contourf(X, Y, nu)
plt.colorbar()
plt.xlabel('x')
plt.ylabel('y')
plt.show()