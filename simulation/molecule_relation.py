"""this module define the relation between 2 molecules"""
from __future__ import annotations
from simulation import DBT1
from electron_coupling.coulomb_matrix import CoulombMatrix
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
 

def create_relation(m1:DBT1.DBT1, m2:DBT1.DBT1,distance:float, Stride:float, translation: tuple[int, int, int])->Relation:
    """This function creates the relation between 2 molecules"""

    CM = gen_cm_matrix(m1, m2, Stride, translation)

    return Relation(distance, CM, translation)

def gen_cm_matrix(m1:DBT1.DBT1, m2:DBT1.DBT1, Stride:float, translation: tuple[int, int, int])->np.ndarray:
    """This function return the coulomb matrix"""
    

    x1 = [a.x for a in m1.atoms]
    y1 = [a.y for a in m1.atoms]
    z1 = [a.z for a in m1.atoms]

    x2 = [a.x + translation[0]*Stride for a in m2.atoms]
    y2 = [a.y + translation[1]*Stride for a in m2.atoms]
    z2 = [a.z + translation[2]*Stride for a in m2.atoms]

    CM_object = CoulombMatrix(x1, x2, y1, y2, z1, z2, DBT1.DBT1.numerator_matrix, DBT1.DBT1.length, False)
    CM = CM_object.gen_CM()

    return CM

