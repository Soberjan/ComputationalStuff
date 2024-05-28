from fenics import *
import matplotlib.pyplot as plt
import os

# msh = Mesh("/home/soberjan/PythonProjects/Fenics/EEG/Debugging/Meshes/donut.xml")
# vtkfile1 = File('msh.pvd')
# vtkfile1 << msh

print(os.environ["DIJITSO_CACHE_DIR"])