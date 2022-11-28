"""This module track down the translation of the repeating space when electron jumping
basic working principle
------
consider a periodic repeating square (when go the the end of left, go reach the right and vice versa)
the square is separated into 4 quadrant, bottom left, bottom right, upper left and upper right

There are 2 things this module can do:
1) define the space when coordinate is given
2) define the possible translation to give the nearest quadrant
"""

class position:
    """store the 3d position of the coordinate"""

    def __init__(self, x: float, y:float, z:float) -> None:
        
        self.x = x
        self.y = y
        self.z = z

class translation:
    """this class track the xyz translation of periodic space"""

    def __init__(self, x: int = 0, y: int = 0, z:int = 0) -> None:
        
        self.x = x
        self.y = y
        self.z = z

def conjugate(translation:translation) -> translation:
    """this method generate the conjugate of the translation object"""

    return translation(-translation.x, -translation.y, -translation.z)

class region_divider:
    """this class define the current quadrant base on the cutting requirement"""

    def __init__(self, x: tuple[float, float], y: tuple[float, float], z: tuple[float, float], stride:int) -> None:
        """arguments:
            x - tuple of x lower bound and upper bound
            y - tuple of y lower bound and upper bound
            z - tuple of z lower bound and upper bound
            stride - steps to get to end of border"""

        self.x = x
        self.y = y
        self.z = z
        self.stride = stride

    def check_boundary(self)->None:
        """this method print the user defined boundary"""

        string_to_print = \
            f"""  |lower|upper\n=================\nx |{self.x[0]}\t|{self.x[1]}\ny |{self.y[0]}\t|{self.y[1]}\nz |{self.z[0]}\t|{self.z[1]}"""

        print(string_to_print)

    def divider(self, axis:str)->float:

        if (axis == 'x'):
            return self.x[1]/self.stride

        elif(axis == 'y'):
            return self.y[1]/self.stride

        else:
            return self.z[1]/self.stride
    
    def get_region(self,coordinate:position)->tuple[int, int, int]:
        """this method find the region of the given coordinate"""

        print(coordinate.x, self.divider('x'), sep = '|')
        x_region = coordinate.x//self.divider('x')
        y_region = coordinate.y//self.divider('y')
        z_region = coordinate.z//self.divider('z')
    
        return int(x_region), int(y_region), int(z_region)


class translation_tracker:
    """This class define the translation space currently in"""
