import numpy as np
import matplotlib.pyplot as plt

# eps = 0.001
# a = 4.0
# b = -0.55
# n = 0.4

eps = 0.02
a = 9.3995
b = -0.405
n = 0.4

def f(Y, t):
    y1, y2 = Y
    return [-y1*(y1-1)*(y1-n)-y2, eps*(y1-a*y2+b)]

xl, xr = -0.25, 1.25
yl, yr = -0.1, 0.2
y1 = np.linspace(xl, xr, 20)
y2 = np.linspace(yl, yr, 20)

Y1, Y2 = np.meshgrid(y1, y2)

t = 0

u, v = np.zeros(Y1.shape), np.zeros(Y2.shape)

NI, NJ = Y1.shape

for i in range(NI):
    for j in range(NJ):
        x = Y1[i, j]
        y = Y2[i, j]
        yprime = f([x, y], t)
        u[i,j] = yprime[0]
        v[i,j] = yprime[1]

Q = plt.quiver(Y1, Y2, u, v, color='r')

plt.xlabel('$y_1$')
plt.ylabel('$y_2$')
plt.xlim([xl, xr])
plt.ylim([yl, yr])
plt.show()