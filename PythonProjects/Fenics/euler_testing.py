import matplotlib.pyplot as plt
import numpy as np


OX = np.array([16, 32, 64, 128, 256, 512, 1024, 2048, 4096])
OY = np.array([0.223231, 0.196017, 0.156569, 0.118207, 0.0864868, 0.0622314, 0.044394, 0.0315309, 0.0223453])
     
plt.plot(OX, OY)
plt.show()
