"""this module define DBT1 molecule that will be used in the project"""

from data_structure.molecule import molecule, atom
import numpy as np
import os
from typing import Optional

class Hydrogen(atom):
    """hydrogen atom, charge = 1"""

    def __init__(self, x: float, y: float, z: float) -> None:
        
        super().__init__(x, y, z, 1)
        self.__name__ = 'Hydrogen'

class Sulphur(atom):
    """Sulphur atom, charge = 16"""

    def __init__(self, x: float, y: float, z: float) -> None:
        super().__init__(x, y, z, 16)
        self.__name__ = 'Sulphur'

class Carbon(atom):
    """Carbon atom, charge = 6"""

    def __init__(self, x: float, y: float, z: float) -> None:
        
        super().__init__(x, y, z, 6)
        self.__name__ = 'Carbon'

class Nitrogen(atom):
    """Nitrogen atom charge = 7"""

    def __init__(self, x: float, y: float, z: float) -> None:

        super().__init__(x, y, z, 7)
        self.__name__ = 'Nitrogen'

class Oxygen(atom):
    """Oxygen atom charge = 8"""

    def __init__(self, x: float, y: float, z: float) -> None:

        super().__init__(x, y, z, 8)
        self.__name__ = 'Oxygen'

def charge(charge1:str,charge2:str)->float:
    """function to convert char to charge value"""
    
    charge = {'H':1,'C':6,'N':7,'O':8,'S':16}
    return charge[charge1]*charge[charge2]

charge_vect = np.vectorize(charge,otypes=[float],cache=True)
"""vectorized charge function"""

def DBT_numerator_matrix(filepath:str)->np.ndarray:
    """
    this method is used to get product of charge 1 and charge 2, as all the molecules are the same,
    there is only one numerator matrix throughout the program
    """
    
    with open(filepath, 'r') as f:

        list_of_elements = f.readlines()
        list_of_elements = [l.rstrip('\n') for l in list_of_elements]

    #this line will generate the matrix
    element_matrix = np.tile(list_of_elements,(len(list_of_elements), 1))
    charge_matrix = charge_vect(element_matrix, element_matrix.T)

    return charge_matrix

class DBT(molecule):
    """DBT1 is the molecule that is used in the simulation, this molecule contains C, H, N, and S"""

    length = 0
    S1 = 0
    S2 = 0
    N = 0
    name = 'DBT'
    molecule_metadata = 'cache/DBT1-molecule-pair'
    reorganization_energy = 0

    def __init__(self, atoms: list[atom], name:str = 'DBT') -> None:
        
        super().__init__(atoms)
        self.__name__ = name

    def __verify__(self) -> None:
        """verify if the molecule constructed is DBT"""

        if len(self.atoms) != self.length :
            print(f'the length of the list is not {self.length}, this is not a {self.name} molecule!')

        if  not isinstance(self.atoms[self.S1],Sulphur):
            print(f'atom {self.S1} is not Sulphur!')

        if not isinstance(self.atoms[self.S2], Sulphur):
            print(f'atom {self.S2} is not Sulphur!')

        if not isinstance(self.atoms[self.N], Nitrogen):
            print(f'atom {self.N} is not Nitrogen!')

    def N_S1_S2_coordinates(self, stride:float, translation:tuple[int, int, int] = (0, 0, 0))-> tuple:
        """This method return the coordinate of the N, S1, and S2 atoms.
        This method is useful when comparing the distance between molecules in periodic space
        this method should only be used together with distance function"""

        x = (self.atoms[self.N].x + translation[0] * stride,
             self.atoms[self.S1].x + translation[0] * stride,
             self.atoms[self.S2].x + translation[0] * stride)

        y = (self.atoms[self.N].y + translation[1] * stride,
             self.atoms[self.S1].y + translation[1] * stride,
             self.atoms[self.S2].y + translation[1] * stride)

        z = (self.atoms[self.N].z + translation[2] * stride,
             self.atoms[self.S1].z + translation[2] * stride,
             self.atoms[self.S2].z + translation[2] * stride)

        return x, y, z
    
    def center_coordinate(self, stride:float, translation:tuple[int, int, int]=(0, 0, 0))->tuple[float,...]:
        
        x = self.atoms[self.N].x + translation[0] * stride
        y = self.atoms[self.N].y + translation[1] * stride
        z = self.atoms[self.N].z + translation[2] * stride

        return x, y, z

def cartesian_distance(x1, x2, y1, y2, z1, z2):
    
    return ((x1-x2)**2 + (y1-y2)**2 + (z1 - z2)**2)**0.5

def DBT_distance(x1:tuple[float, float, float],
                  y1:tuple[float, float, float],
                  z1:tuple[float, float, float],
                  x2:tuple[float, float, float],
                  y2:tuple[float, float, float],
                  z2:tuple[float, float, float])->float:

    """This method calculates the distance between two DBT molecule"""
    
    distance_N_N = cartesian_distance(x1[0], x2[0], y1[0], y2[0], z1[0], z2[0])
    distance_N_S1 = cartesian_distance(x1[0], x2[1], y1[0], y2[1], z1[0], z2[1])
    distance_N_S2 = cartesian_distance(x1[0], x2[2], y1[0], y2[2], z1[0], z2[2])


    distance_S1_N = cartesian_distance(x1[1], x2[0], y1[1], y2[0], z1[1], z2[0])
    distance_S1_S1 = cartesian_distance(x1[1], x2[1], y1[1], y2[1], z1[1], z2[1])
    distance_S1_S2 = cartesian_distance(x1[1], x2[2], y1[1], y2[2], z1[1], z2[2])

    distance_S2_N = cartesian_distance(x1[2], x2[0], y1[2], y2[0], z1[2], z2[0])
    distance_S2_S1 = cartesian_distance(x1[2], x2[1], y1[2], y2[1], z1[2], z2[1])
    distance_S2_S2 = cartesian_distance(x1[2], x2[2], y1[2], y2[2], z1[2], z2[2])


    Cut_off_distance =(distance_N_S1 + distance_N_S2 + distance_N_N + \
                        distance_S1_S1 + distance_S1_S2 + distance_S1_N + \
                        distance_S2_S1 + distance_S2_S2 + distance_S2_N)/9

    return Cut_off_distance

