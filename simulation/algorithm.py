"""This module contains the algorithm that is used by this project only"""
from typing import Optional, Union
from simulation import DBT1
from simulation.molecule_relation import create_relation, Relation
from itertools import product
from functools import lru_cache

def molecule_is_cut(m:DBT1.DBT1, limit)->tuple[bool, bool, bool]:
    """this function detects if the molecule is cut by boundary, use for this project only"""

    x = [a.x for a in m.atoms]
    y = [a.y for a in m.atoms]
    z = [a.z for a in m.atoms]
    
    x_is_cut, y_is_cut, z_is_cut = False, False, False

    if (max(x) - min(x) > limit):
        
        x_is_cut = True
    
    if (max(y) - min(y) > limit):

        y_is_cut = True

    if (max(z) - min(z) > limit):

        z_is_cut = True

    return x_is_cut, y_is_cut, z_is_cut

def complete_molecule(m:DBT1.DBT1, boundary)->DBT1.DBT1:
    """This function make the molecule complete by fixing the x, y, z coordinate of the atoms that is cut by border"""
    
    is_cut = molecule_is_cut(m, 0.5 * boundary)
    if any(is_cut):

        if is_cut[0]:

            for index, atom in enumerate(m.atoms):

                if atom.x < 5:

                    m.atoms[index].x += boundary

        if is_cut[1]:

            for index, atom in enumerate(m.atoms):

                if atom.y < 5:

                    m.atoms[index].y += boundary

        if is_cut[2]:

            for index, atom in enumerate(m.atoms):

                if atom.z < 5:

                    m.atoms[index].z += boundary

        return DBT1.DBT1(m.atoms, m.__name__)

    return m

@lru_cache(maxsize=1)
def distance(m1:DBT1.DBT1, m2:DBT1.DBT1, Stride, Translation:tuple[int, int, int])->float:

    Coord1, Coord2 = m1.N_S1_S2_coordinates(Stride), m2.N_S1_S2_coordinates(stride = Stride, translation = Translation)

    return DBT1.DBT1_distance(*Coord1, *Coord2)

def DBT1_pair_is_close(m1:DBT1.DBT1, m2: DBT1.DBT1, cut_off:float, Stride, Translation:tuple[int, int, int] = (0, 0, 0))->bool:
    """This function tell whether the two molecule is close enough to be considered as pair
        
        return: bool"""
    
    molecular_distance = distance(m1, m2, Stride, Translation)

    if (molecular_distance<= cut_off):

        return True

    return False

def _translation_checker(x1:tuple[float, float, float], y1:tuple[float, float, float], z1:tuple[float, float, float],
                                 x2:tuple[float, float, float], y2:tuple[float, float, float], z2:tuple[float, float, float],
                                 stride:float)->tuple[int, int, int]:
    """This function should return the most possible transformation for the nearest distance in the periodic space"""

    #check the 9 vectors of x, then check the 9 vectors of y, then z
    
    def vec_checker(tup1:tuple[float, float, float], tup2:tuple[float, float, float]):
        """
        this function check if there is a possible transformation that can reduce the distance in the periodic space,
        if there is no way to reduce anymore, then 0 will be return, 
        if a negative transformation is possible then return -1, and vice versa
        """
        avg_dist = 0

        for i in product(tup2, tup1):
            
            vect = i[0] - i[1]

            if -0.5*stride <= vect <= 0.5*stride: #if any of the vector fall into this region the check can be stopped immediately, return no transformation
            
                return 0
            
            avg_dist += vect/9

        if avg_dist > 0:

            return -1

        return 1
    
    x_transformation = vec_checker(x1, x2)
    y_transformation = vec_checker(y1, y2)
    z_transformation = vec_checker(z1, z2)

    return x_transformation, y_transformation, z_transformation

def periodic_translation(m1:DBT1.DBT1, m2:DBT1.DBT1, Stride)->tuple[int, int, int]:
    """This function is a helper/wrapper function to pass arguments into the transformation function"""

    Coord1, Coord2 = m1.N_S1_S2_coordinates(Stride), m2.N_S1_S2_coordinates(Stride)
    return _translation_checker(*Coord1, *Coord2, Stride)

def molecular_pair_relation(frm:DBT1.DBT1, to:DBT1.DBT1, cut_off:float, Stride:float)->Union[Relation, None]:
    """This function should tell the relationship between 2 molecules,
    if there is no relation then it should return None"""

    translation = periodic_translation(frm, to, Stride)
    if DBT1_pair_is_close(frm, to, cut_off, Stride, translation):
        
        Distance = distance(frm, to, Stride, translation)
        return create_relation(frm, to, Distance, Stride, translation)
    
    return None

def molecular_pair_relation_brute_force_check(frm:DBT1.DBT1, to:DBT1.DBT1, cut_off:float, Stride:float)->Union[Relation, None]:
    """This function should tell the relationship between 2 molecules,
    if there is no relation then it should return None"""

    translation = periodic_translation(frm, to, Stride)

    if DBT1_pair_is_close(frm, to, cut_off, Stride, (0, 0, 0)):

        Distance = distance(frm, to, Stride, translation)
        return create_relation(frm, to, Distance, Stride,(0, 0, 0))
    
    #return None
    
    tup = (-1, 0, 1)
    for i in product(tup, tup, tup):

        if DBT1_pair_is_close(frm, to, cut_off, Stride, i):

            Distance = distance(frm, to, Stride, i)
            return create_relation(frm, to, Stride,Distance, i)
    
    return None
