"""This module is the simulation of electron transfer"""
from frame_data_processing.frame import frame

class simulation:
    """This class define how the simulation should be run from a single frame"""

    def __init__(self, frame:frame) -> None:
        """initialization, take a frame as argument so that all the data can be access from the simulation object
        for the simulation, stopping condition should also be define"""
        
        self.Frame = frame


