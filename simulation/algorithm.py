"""This module contains the algorithm that is used by this project only"""
from typing import Union
from simulation import DBT1
from simulation.molecule_relation import Relation
from itertools import product

def molecule_is_cut(m:DBT1.DBT1)->tuple[bool, bool, bool]:
    """this function detects if the molecule is cut by boundary, use for this project only"""

    limit = 5 #5 nm
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

def complete_molecule(m:DBT1.DBT1)->DBT1.DBT1:
    """This function make the molecule complete by fixing the x, y, z coordinate of the atoms that is cut by border"""
    
    boundary = 9.66968

    is_cut = molecule_is_cut(m)
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

        return DBT1.DBT1(m.atoms)

    return m

def DBT1_pair_is_close(m1:DBT1.DBT1, m2: DBT1.DBT1, cut_off:float, Stride, Translation:tuple[int, int, int] = (0, 0, 0))->bool:
    """This function tell whether the two molecule is close enough to be considered as pair
        
        return: bool"""
    
    Coord1, Coord2 = m1.N_S1_S2_coordinates(Stride), m2.N_S1_S2_coordinates(stride = Stride, translation = Translation)

    molecular_distance = DBT1.DBT1_distance(*Coord1, *Coord2)

    if (molecular_distance<= cut_off):

        return True

    return False

def molecular_pair_relation(frm:DBT1.DBT1, to:DBT1.DBT1, cut_off:float, Stride:float)->Union[Relation, None]:
    """This function should tell the relationship between 2 molecules,
    if there is no relation then it should return None"""

    if DBT1_pair_is_close(frm, to, cut_off, Stride):

        return Relation(frm, to, (0, 0, 0))

    return None

    #this part will check the transformation
    tup = (-1, 0, 1)
    ls = [i for i in product(tup, tup, tup) if i not in [(0, 0, 0)]]

    for i in ls:

        if DBT1_pair_is_close(frm, to, cut_off, Stride, i):

            print(f"pair found at translation {i}")
            return Relation(frm, to, i)

    return None
