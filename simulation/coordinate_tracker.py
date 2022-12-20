"""
This module contains coordinate tracker class to keep track of the trajectory of electron transport
"""

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

def Plot_trajectory(x:np.ndarray, y:np.ndarray, z:np.ndarray):

    fig = plt.figure()

    ax = fig.add_subplot(111, projection='3d')

    ax.plot(x,y,z)

    fig.savefig('testing.png')

