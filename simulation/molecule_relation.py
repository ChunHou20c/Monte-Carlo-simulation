"""this module define the relation between 2 molecules"""
from __future__ import annotations
from simulation import DBT1
from electron_coupling.pair_coulomb_matrix import pair_coulomb_matrix
from electron_coupling import coulomb_matrix
import numpy as np

class Relation:
    """This class stores the relationship between 2 molecule"""

    def __init__(self, distance:float, Coulomb_matrix:np.ndarray, translation:tuple[int, int, int] = (0, 0, 0)) -> None:
        
        self.distance = distance
        self.coulomb_matrix = Coulomb_matrix
        self.translation = translation


    def conjugate(self)->Relation:
        """This method return the conjugate of the relation"""

        translation = tuple(-i for i in self.translation)
        return Relation(self.distance, self.coulomb_matrix, translation)
 

def create_relation(m1:DBT.DBT, m2:DBT.DBT,distance:float, Stride:float, translation: tuple[int, int, int])->Relation:
    """This function creates the relation between 2 molecules"""

    CM = gen_cm_matrix(m1, m2, Stride, translation)

    return Relation(distance, CM, translation)

def gen_cm_matrix(m1:DBT.DBT, m2:DBT.DBT, Stride:float, translation: tuple[int, int, int])->np.ndarray:
    """This function return the coulomb matrix"""
    

    x1 = [a.x for a in m1.atoms]
    y1 = [a.y for a in m1.atoms]
    z1 = [a.z for a in m1.atoms]

    x2 = [a.x + translation[0]*Stride for a in m2.atoms]
    y2 = [a.y + translation[1]*Stride for a in m2.atoms]
    z2 = [a.z + translation[2]*Stride for a in m2.atoms]

    x = x1 + x2
    y = y1 + y2
    z = z1 + z2

    numerator_matrix = coulomb_matrix.gen_numerator_matrix(m1.molecule_metadata)
    #print(numerator_matrix.shape)
    #print(type(m1))
    #print(type(m2))

    CM_object = pair_coulomb_matrix(x, y, z, numerator_matrix)
    CM = CM_object.gen_full_coulomb_matrix()

    return CM
