import matplotlib.pyplot as plt
import numpy as np

def f(x):
    if x >= 2:
        return f(x-2)
    if 0 <= x <= 0.33333:
        return 0
    if 0.33 < x <= 0.666666:
        return 3 * x - 1
    if 0.66 < x <= 1:
        return 1
    if 1 < x <= 2:
        return -x + 2

OX = np.linspace(0, 2, 10000)

# OY = []
# for x in OX:
#     OY.append(f(x))
# plt.plot(OX, OY)

OY1 = []
for x in OX:
    OY1.append(0.5 * f(3 * x))
#plt.plot(OX, OY1)

OY2 = []
for x in OX:
    OY2.append(0.25 * f(27 * x))
#plt.plot(OX, OY2)

OY3 = []
for x in OX:
    OY3.append(0.125 * f(81 * x))


OY4 = []
for x in OX:
    OY4.append(1/16 * f(81 * x))
    
OYn = []
for x in range(len(OY1)):
    OYn.append(OY1[x] + OY2[x] + OY3[x] + OY4[x])
plt.plot(OX, OYn)

plt.show()