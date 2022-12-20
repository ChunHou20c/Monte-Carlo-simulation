"""
This module contain the multiframe (dynamic) simulation

the dynamic simulation is just a class that stores multiple simulation object and provide methods to traverse
through the simulation list
"""

from simulation import simulation
from simulation import algorithm
from simulation import coordinate_tracker
import os
import random
import numpy as np

class dynamic_simulation:

    def __init__(self, static_frames:list[simulation.Simulation]) -> None:
        
        self.frames_dict:dict[float, simulation.Simulation] = build_dict(static_frames)
        self.timeframe = 0
        self.total_jump = 10000
        self.coordinates = [[],[],[]]
        self.vector = np.array((0.0000,0.0000,0.0000))

    def inspect(self):

        print(self.frames_dict.keys())

    def run(self):
        """This is the main loop for the simulation"""

        current_frame = self.frames_dict[self.timeframe]
        keys = list(current_frame.graph.get_vertices())
        current_key = random.choice(keys)
        
        initial_x, initial_y, initial_z = current_frame.graph.get_vertex(current_key).molecule.center_coordinate(current_frame.box_width, (0,0,0))
        
        self.coordinates[0].append(initial_x)
        self.coordinates[1].append(initial_y)
        self.coordinates[2].append(initial_z)

        for _ in range(self.total_jump):
            
            next_frame = algorithm.closest_value(list(self.frames_dict.keys()), self.timeframe)
            current_frame = self.frames_dict[next_frame]
            
            current_key, jumping_time, vector = current_frame.single_jump(current_key)
            
            current_x, current_y, current_z = np.array([self.coordinates[0][-1], self.coordinates[1][-1], self.coordinates[2][-1]]) + np.array(vector)
            
            self.coordinates[0].append(current_x)
            self.coordinates[1].append(current_y)
            self.coordinates[2].append(current_z)
            
            self.timeframe += jumping_time*1e12
            print(f'{self.timeframe = }')
            self.vector += np.array(vector)

    def export_trajectory(self):
        """This method export the trajectory data in numpy array to draw the trajectory"""

        return np.array(self.coordinates[0]), np.array(self.coordinates[1]), np.array(self.coordinates[2])

    def Plot_trajectory(self):
        """This method plot the trajectory and save it into a file"""

        coordinate_tracker.Plot_trajectory(*self.export_trajectory())

def build_dict(static_frames:list[simulation.Simulation])->dict:
    
    dict_to_return = {}

    for frame in static_frames:

        dict_to_return[frame.timestamp] = frame

    return dict_to_return

def build_dynamic_simulation(dir:str):
    """This function build the dynamic simulation object from a directory"""

    sim_list = []
    for file in os.listdir(dir):
        full_path = os.path.join(dir, file)
        sim = simulation.Simulation(full_path)
        sim_list.append(sim)

    Dynamic_sim = dynamic_simulation(sim_list)

    return Dynamic_sim
