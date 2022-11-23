"""this module define the relation between 2 molecules"""
from __future__ import annotations
from simulation import DBT1


class Relation:
    """This class stores the relationship between 2 molecule"""

    def __init__(self, m1:DBT1.DBT1, m2:DBT1.DBT1, translation:tuple[int, int, int] = (0, 0, 0)) -> None:
        
        self.m1 = m1
        self.m2 = m2
        self.translation = translation

    def conjugate(self)->Relation:
        """This method return the conjugate of the relation"""

        translation = tuple(-i for i in self.translation)
        return Relation(self.m2, self.m1, translation)
        
