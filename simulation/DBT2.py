"""this module define DBT1 molecule that will be used in the project"""

from data_structure.molecule import molecule, atom
import numpy as np
import os
from typing import Optional

from simulation import DBT
#64 atoms

class DBT2(DBT.DBT):
    """DBT1 is the molecule that is used in the simulation, this molecule contains 64 atoms and formed from C, H, N, and S"""

    length = 64
    S1 = 27
    S2 = 5
    N = 20
    molecule_metadata = 'cache/DBT2-molecule-pair'
    reorganization_energy = 0.373

    def __init__(self, atoms: list[atom], name:str = 'DBT1') -> None:
        
        super().__init__(atoms)
        self.__name__ = name
