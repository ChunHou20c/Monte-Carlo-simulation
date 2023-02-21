"""
This module contains coordinate tracker class to keep track of the trajectory of electron transport
"""

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import imageio
import seaborn as sns
import pandas as pd

def Plot_trajectory(x:np.ndarray, y:np.ndarray, z:np.ndarray, t:float, name:str):

    fig = plt.figure(figsize=(6,6))
    
    ax = fig.add_subplot(111, projection='3d')

    ax.plot(x,y,z, color='grey')
    ax.plot(x[-1],y[-1], z[-1], color = 'black', marker='o')
    ax.set_zlabel('z(nm)', fontsize=14)
    ax.set_xlabel('x (nm)', fontsize=14)
    ax.set_ylabel('y (nm)', fontsize=14)
    plt.suptitle('Electron Trajectory', fontsize=14)
    plt.title(f'step = {len(x)}, t = {t:.3f} ps', fontsize = 12)
    fig.savefig(name)
    
    plt.close(fig)

def create_trajectory_frames(x:np.ndarray, y:np.ndarray, z:np.ndarray, t:np.ndarray, dir:str):
    """
    This function create the frames for gif generation
    """

    for i, _ in enumerate(x, 1):

        frame_name = f'{dir}/{dir}_{i}.png'

        Plot_trajectory(x[:i], y[:i], z[:i],t[i-1],frame_name)

    frames = []
    for i, _ in enumerate(x,1):

        image = imageio.v2.imread(f'{dir}/{dir}_{i}.png')
        frames.append(image)

    imageio.mimsave('electron_trajectory_10fps.gif', frames, fps = 10)
    imageio.mimsave('electron_trajectory_25fps.gif', frames, fps = 25)

def plot_square_distance(x:np.ndarray, t:np.ndarray, name:str):
    """
    This function plot the square distance vs time
    """
    
    #make the data into dataframe

    sd_column = 'square distance (nm^2)'
    time_coloumn = 'time (ps)'

    df = pd.DataFrame(list(zip(x, t)), columns=[sd_column, time_coloumn])
    sns.set(rc={'figure.figsize':(6, 6)})
    sns.set_style('dark')
    plot = sns.lineplot(x=time_coloumn, y=sd_column, data=df, color='blue')
    plot.axes.set_title('Square Displacement vs Time', fontsize=14)
    plot.set_xlabel("time (ps)",fontsize=13)
    plot.set_ylabel("squared distance ($nm^2$)",fontsize=13)
    fig = plot.get_figure()
    fig.savefig(name)
    
    plt.close(fig)

def create_msd_gif(x:np.ndarray, t:np.ndarray, dir:str):
    """
    This function create gif of msd
    """

    for i, _ in enumerate(x, 1):

        frame_name = f'{dir}/{dir}_{i}.png'

        plot_square_distance(x[:i], t[:i],frame_name)

    frames = []
    for i, _ in enumerate(x,1):

        image = imageio.v2.imread(f'{dir}/{dir}_{i}.png')
        frames.append(image)

    imageio.mimsave('msd_10fps.gif', frames, fps = 10)
    imageio.mimsave('msd_25fps.gif', frames, fps = 25)
