import math
import numpy as np
import matplotlib.pyplot as plt

L = 15
dt = 5
Du = 5
def U(x, t):
    return np.cos(np.pi*x/L) * 5*np.exp(-Du*np.pi**2*t/L**2)

x = np.linspace(0, L, 100)
t = 0
while t <= 10:
    u = U(x, t)
    plt.plot(x, u, label=t)
    t += dt

plt.xlabel('x')
plt.ylabel('u')
plt.legend()
plt.show()