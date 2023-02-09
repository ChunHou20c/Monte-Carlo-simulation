"""this module define DBT1 molecule that will be used in the project"""

from data_structure.molecule import molecule, atom
from simulation import DBT
import numpy as np
import os
from typing import Optional

class DBT1(DBT.DBT):
    """DBT1 is the molecule that is used in the simulation, this molecule contains 56 atoms and formed from C, H, N, and S"""

    length = 56
    S1 = 42
    S2 = 5
    N = 20
    name = 'DBT1'
    molecule_metadata = 'cache/DBT1-molecule-pair'
    reorganization_energy = 0.180

    def __init__(self, atoms: list[atom], name:str = 'DBT1') -> None:
        
        super().__init__(atoms)
        self.__name__ = name

