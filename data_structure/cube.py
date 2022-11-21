"""This module define a cube in cartesian coordinate"""
from __future__ import annotations

class Cube:
    """A cube object have a lower boundary and side length"""

    def __init__(self, side:float, x_lower_bound:float = 0, y_lower_bound:float = 0, z_lower_bound:float = 0) -> None:
        
        self.x_lower_bound = x_lower_bound
        self.y_lower_bound = y_lower_bound
        self.z_lower_bound = z_lower_bound
        self.side = side

    def __str__(self) -> str:
        
        x_range = (self.x_lower_bound, self.x_lower_bound + self.side)
        y_range = (self.y_lower_bound, self.y_lower_bound + self.side)
        z_range = (self.z_lower_bound, self.z_lower_bound + self.side)
        
        return f"{x_range=}, {y_range=}, {z_range=}"

    def neighbor_cube(self, relative_x:int=0, relative_y:int=0, relative_z:int=0)->Cube:
        """This method return the neighboring cube with is relative to this cube"""

        return Cube(self.side, 
                    self.x_lower_bound + relative_x * self.side,
                    self.y_lower_bound + relative_y * self.side,
                    self.z_lower_bound + relative_z * self.side)
        

        

