import math
import numpy as np
import matplotlib.pyplot as plt

L = 15
dt = 5
Du = 5.0
eps = 0.01
a = 1.0
B = Du*np.pi**2/L**2
A = 5.0/(eps * a - B)
def V(x, t):
    return eps*A*np.cos(np.pi*x/L)*(np.exp(-B*t)-np.exp(-eps*a*t))

x = np.linspace(0, L, 100)
t = 0
while t <= 10:
    v = V(x, t)
    plt.plot(x, v, label=t)
    t += dt

plt.xlabel('x')
plt.ylabel('v')
plt.legend()
plt.show()