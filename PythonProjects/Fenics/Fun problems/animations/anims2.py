import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import PillowWriter

fig, axis = plt.subplots(2, 1)
l1, = axis[0].plot([], [], 'k-')
l2, = axis[1].plot([], [], 'k-')
axis[0].set_xlim(-5, 5)
axis[0].set_ylim(-5, 5)
axis[1].set_xlim(-5, 5)
axis[1].set_ylim(-5, 5)

def func(x):
    return np.sin(x)*3
def func2(x):
    return np.cos(x)*3

writer = PillowWriter(fps=15)

xlist = np.linspace(-5, 5, 100)
ylist = func(xlist)
ylist2 = func2(xlist)

with writer.saving(fig, 'SineWave.gif', 100):
    for i in range(1, 10):
        l1.set_data(xlist, ylist)
        l2.set_data(xlist, ylist2)
        writer.grab_frame()
        ylist = func(xlist * i)
        ylist2 = func2(xlist * i)
