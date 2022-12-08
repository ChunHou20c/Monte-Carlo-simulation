"""
This module contain the multiframe (dynamic) simulation

the dynamic simulation is just a class that stores multiple simulation object and provide methods to traverse
through the simulation list
"""

from simulation import simulation
import os

class dynamic_simulation:

    def __init__(self, static_frames:list[simulation.Simulation]) -> None:
        
        self.frames_dict = build_dict(static_frames)
    
def build_dict(static_frames:list[simulation.Simulation])->dict:
    
    dict_to_return = {}

    for frame in static_frames:

        dict_to_return[frame.timestamp] = frame

    return dict_to_return

def build_dynamic_simulation(dir:str):
    """This function build the dynamic simulation object from a directory"""

    sim_list = []
    for file in os.listdir(dir):

        sim = simulation.Simulation(f'{dir}/{file}')
        sim_list.append(sim)

    Dynamic_sim = dynamic_simulation(sim_list)

    return Dynamic_sim
